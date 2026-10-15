# Pass 1: Select a near-minimal set of templates, such that the depth of the selected
# templates is at least min(D, d(p)) at every position p, where d(p) is the full depth.
# See algorithm.md for the algorithm.

# Orders the blocks in `pending` and `active` by their first field alone (start, or stop).
# Comparing whole tuples would also compare the template on ties, and check equality of
# the first field before `<`. Ties don't affect the selection, so the order of equal
# keys doesn't matter.
const BY_FIRST = By(first)

"""
    template_hash(name::Union{AbstractString, Nothing})::UInt64

Hash a QNAME, where `nothing` is the missing QNAME `*`. Templates are identified by
this hash. A collision merges two templates into one selection unit, which costs
optimality, but not correctness.

The `StringView` returned by `query_name(::BAMRecord)` hashes its bytes in place with
Base's string hash, so this doesn't allocate, and equals the hash of the same `String`.
"""
template_hash(name::AbstractString)::UInt64 = hash(name) % UInt64
template_hash(::Nothing)::UInt64 = template_hash("*")

"""
    Sweep(depth::Int)

State of the pass-1 sweep over the blocks of one reference at a time, selecting templates
such that the selected depth is at least `min(depth, d(p))` at every position `p`.

Positions are only evaluated at events (block starts, and one past block ends),
since the depths are constant between them.
"""
mutable struct Sweep
    const depth::Int
    # Hashes of selected templates. Persists across references.
    const selected::Set{UInt64}
    # Blocks not yet reached by the sweep: (start, stop, template)
    const pending::Vector{Tuple{Int, Int, UInt64}}
    # Blocks overlapping the current position, as a bucket queue keyed by stop: a ring of
    # linked lists of entries, where the blocks ending at `stop` are in bucket
    # `(stop & (length(heads) - 1)) + 1`. All active stops are in
    # `base:(base + length(heads) - 1)`, so each bucket holds blocks of a single stop.
    # The ring grows when a longer block is activated.
    heads::Vector{Int32}
    # The next entry of each entry's list, or 0 at the end of a list.
    # For free entries, the next free entry.
    const entry_next::Vector{Int32}
    const entry_template::Vector{UInt64}
    # First free entry, or 0 if none
    free::Int32
    # The position of the last event, or less
    base::Int
    # A lower bound of the stops of the active blocks, which is at least `base`
    lo::Int
    # Number of active blocks, which is the full depth
    n_active::Int
    # Active blocks of unselected templates, as a max-heap: (stop, template).
    # Entries are removed lazily, so it may contain blocks that have ended, or whose
    # template has been selected.
    const candidates::Vector{Tuple{Int, UInt64}}
    # For each unselected template with active blocks, its number of active blocks.
    # Hence, an active block's template is selected iff it's not a key here.
    const unselected::Dict{UInt64, Int}
    # Number of active blocks belonging to unselected templates
    n_unselected::Int
    # Number of active blocks belonging to selected templates
    selected_depth::Int
end

function Sweep(depth::Int)
    depth > 0 || throw(ArgumentError("depth must be positive"))
    return Sweep(
        depth, Set{UInt64}(), Tuple{Int, Int, UInt64}[], zeros(Int32, 1024), Int32[],
        UInt64[], Int32(0), 0, 0, 0, Tuple{Int, UInt64}[], Dict{UInt64, Int}(), 0, 0,
    )
end

function add_block!(sweep::Sweep, block::UnitRange{Int}, template::UInt64)
    heappush!(sweep.pending, (first(block), last(block), template), BY_FIRST)
    return sweep
end

# The smallest stop of the active blocks, if any. Raises `lo` to it.
function min_active_stop!(sweep::Sweep)::Union{Int, Nothing}
    sweep.n_active == 0 && return nothing
    heads = sweep.heads
    mask = length(heads) - 1
    lo = sweep.lo
    while (@inbounds heads[(lo & mask) + 1]) == 0
        lo += 1
    end
    sweep.lo = lo
    return lo
end

# The next position where a block starts or ends, if any
function next_event(sweep::Sweep)::Union{Int, Nothing}
    pending = sweep.pending
    stop = min_active_stop!(sweep)
    if isempty(pending)
        stop === nothing ? nothing : stop + 1
    elseif stop === nothing
        (@inbounds pending[1])[1]
    else
        min((@inbounds pending[1])[1], stop + 1)
    end
end

"""
    sweep_until!(sweep::Sweep, limit::Int)

Process all events at positions before `limit`. The caller must not add blocks that
start before `limit` afterwards.
"""
function sweep_until!(sweep::Sweep, limit::Int)
    while true
        position = next_event(sweep)
        (position === nothing || position >= limit) && break
        process_event!(sweep, position)
    end
    return sweep
end

"""
    finish_reference!(sweep::Sweep)

Process all remaining events. Afterwards, `sweep` is ready for the next reference.
All blocks have ended, so the entries left in `candidates` are stale, and are discarded.
"""
function finish_reference!(sweep::Sweep)
    sweep_until!(sweep, typemax(Int))
    empty!(sweep.candidates)
    # Positions restart on the next reference
    sweep.base = 0
    sweep.lo = 0
    return sweep
end

# Grow the ring of `active` so it can hold a block ending at `stop`
@noinline function grow_ring!(sweep::Sweep, stop::Int)
    old = sweep.heads
    mask = length(old) - 1
    new = zeros(Int32, nextpow(2, 2 * (stop - sweep.base + 1)))
    new_mask = length(new) - 1
    for s in sweep.base:(sweep.base + mask)
        @inbounds new[(s & new_mask) + 1] = old[(s & mask) + 1]
    end
    sweep.heads = new
    return nothing
end

function process_event!(sweep::Sweep, position::Int)
    (; pending, unselected, candidates, entry_next, entry_template) = sweep
    # Remove blocks that ended before this position
    heads = sweep.heads
    mask = length(heads) - 1
    for stop in sweep.lo:(position - 1)
        sweep.n_active == 0 && break
        i = (stop & mask) + 1
        entry = @inbounds heads[i]
        entry == 0 && continue
        @inbounds heads[i] = 0
        while entry != 0
            template = @inbounds entry_template[entry]
            n = get(unselected, template, 0)
            if n == 0
                sweep.selected_depth -= 1
            else
                sweep.n_unselected -= 1
                if n == 1
                    delete!(unselected, template)
                else
                    unselected[template] = n - 1
                end
            end
            sweep.n_active -= 1
            next = @inbounds entry_next[entry]
            @inbounds entry_next[entry] = sweep.free
            sweep.free = entry
            entry = next
        end
    end
    sweep.lo = max(sweep.lo, position)
    sweep.base = position
    # Activate blocks starting at this position
    while !isempty(pending) && (@inbounds pending[1])[1] <= position
        (_, stop, template) = heappop!(pending, BY_FIRST)
        stop >= position + length(sweep.heads) && grow_ring!(sweep, stop)
        entry = sweep.free
        if entry == 0
            push!(entry_next, 0)
            push!(entry_template, template)
            entry = Int32(length(entry_next))
        else
            sweep.free = @inbounds entry_next[entry]
            @inbounds entry_template[entry] = template
        end
        heads = sweep.heads
        i = (stop & (length(heads) - 1)) + 1
        @inbounds entry_next[entry] = heads[i]
        @inbounds heads[i] = entry
        sweep.n_active += 1
        sweep.lo = min(sweep.lo, stop)
        if template in sweep.selected
            sweep.selected_depth += 1
        else
            unselected[template] = get(unselected, template, 0) + 1
            sweep.n_unselected += 1
            heappush!(candidates, (stop, template), Reverse)
        end
    end
    # Select templates until the target depth is reached, preferring the template
    # whose block extends the furthest.
    # Every active block of an unselected template has an entry in `candidates`,
    # so it can't run empty while selected_depth < n_active.
    target = min(sweep.depth, sweep.n_active)
    while sweep.selected_depth < target
        (stop, template) = heappop!(candidates, Reverse)
        stop < position && continue
        n = pop!(unselected, template, 0)
        n == 0 && continue
        push!(sweep.selected, template)
        sweep.selected_depth += n
        sweep.n_unselected -= n
    end
    # Remove stale entries, so `candidates` doesn't grow without bound in regions
    # where no selection is needed. This doesn't use `filter!`, which would shrink the
    # vector's memory, so it would be reallocated as it grows again.
    if length(candidates) > 2 * sweep.n_unselected + 1024
        j = 0
        for candidate in candidates
            (stop, template) = candidate
            if stop >= position && haskey(unselected, template)
                j += 1
                @inbounds candidates[j] = candidate
            end
        end
        resize!(candidates, j)
        heapify!(candidates, Reverse)
    end
    return sweep
end

# TODO: More efficient to return a BitVector? Must be stored as set of hashes since
# reads may appear several times and must be deduplicated, but we can map these to
# their first occurrence in the file and then store a BitVector of kept reads.
# This makes lookup in second pass faster, and takes less memory in second pass
# (though more in first pass)
"""
    select_templates(reader::BAMReader, depth::Int)::Set{UInt64}

Read all remaining records of `reader`, and return the [`template_hash`](@ref)es of a
near-minimal set of templates, such that at every reference position `p`, the depth of
the selected templates is at least `min(depth, d(p))`, where `d(p)` is the depth of
all mapped records.
Throw an `OverburdenError` if the records are not coordinate sorted, or have a
malformed CIGAR, or an `XAMError` if the file is otherwise malformed.
"""
function select_templates(reader::BAMReader, depth::Int)::Set{UInt64}
    sweep = Sweep(depth)
    checker = SortChecker(length(XAM2.references(reader)))
    # Index of the reference whose blocks are in `sweep`, or 0 if none
    sweep_reference = 0
    n = 0
    while true
        record = XAM2.try_read!(reader)
        record === nothing && break
        record isa XAMError && throw(record)
        n += 1
        check_sorted!(checker, record, n)
        is_mapped(record) || continue
        index = @something reference_index(record) continue
        position = @something left_position(record) continue
        if index != sweep_reference
            finish_reference!(sweep)
            sweep_reference = index
        else
            sweep_until!(sweep, position)
        end
        cig = cigar(record)
        cig === nothing && continue
        cig isa BAMCIGAR || throw_bad_cigar(record, n, cig)
        template = template_hash(query_name(record))
        # Blocks cover M, =, X and D operations, and are split by N.
        # They are relative to the record's position.
        for block in aligned_blocks(CIGARStrings.ref, CIGARStrings.aln, cig)
            add_block!(sweep, block .+ (position - 1), template)
        end
    end
    finish_reference!(sweep)
    return sweep.selected
end

function select_templates(path::String, depth::Int; n_workers::Int = Threads.nthreads())::Set{UInt64}
    reader = open_bam(path; n_workers, copy = false)
    return try
        select_templates(reader, depth)
    finally
        close(reader)
    end
end

@noinline function throw_bad_cigar(record::BAMRecord, n::Int, err::Union{CIGARError, XAMError})
    name = something(query_name(record), "*")
    throw(OverburdenError("malformed CIGAR in record $n (\"$name\"): " * sprint(showerror, err)))
end
