# Pass 1: Select a near-minimal set of templates, such that the depth of the selected
# templates is at least min(D, d(p)) at every position p, where d(p) is the full depth.
# See algorithm.md for the algorithm.

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
    # Blocks overlapping the current position: (stop, template).
    # Its length is the full depth.
    const active::Vector{Tuple{Int, UInt64}}
    # Active blocks of unselected templates, as a max-heap: (stop, template).
    # Entries are removed lazily, so it may contain blocks that have ended, or whose
    # template has been selected.
    const candidates::Vector{Tuple{Int, UInt64}}
    # For each unselected template with active blocks, the stops of its active blocks.
    # Hence, an active block's template is selected iff it's not a key here.
    const unselected::Dict{UInt64, Vector{Int}}
    # Emptied vectors from `unselected`, for reuse
    const spare::Vector{Vector{Int}}
    # Number of active blocks belonging to unselected templates
    n_unselected::Int
    # Number of active blocks belonging to selected templates
    selected_depth::Int
end

function Sweep(depth::Int)
    depth > 0 || throw(ArgumentError("depth must be positive"))
    return Sweep(
        depth, Set{UInt64}(), Tuple{Int, Int, UInt64}[], Tuple{Int, UInt64}[],
        Tuple{Int, UInt64}[], Dict{UInt64, Vector{Int}}(), Vector{Int}[], 0, 0,
    )
end

function add_block!(sweep::Sweep, block::UnitRange{Int}, template::UInt64)
    heappush!(sweep.pending, (first(block), last(block), template))
    return sweep
end

# The next position where a block starts or ends, if any
function next_event(sweep::Sweep)::Union{Int, Nothing}
    (; pending, active) = sweep
    if isempty(pending)
        isempty(active) ? nothing : (@inbounds active[1])[1] + 1
    elseif isempty(active)
        (@inbounds pending[1])[1]
    else
        min((@inbounds pending[1])[1], (@inbounds active[1])[1] + 1)
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
    return sweep
end

function process_event!(sweep::Sweep, position::Int)
    (; active, pending, unselected, candidates) = sweep
    # Remove blocks that ended before this position
    while !isempty(active) && (@inbounds active[1])[1] < position
        (stop, template) = heappop!(active)
        stops = get(unselected, template, nothing)
        if stops === nothing
            sweep.selected_depth -= 1
        else
            i = findfirst(==(stop), stops)::Int
            @inbounds stops[i] = stops[end]
            pop!(stops)
            sweep.n_unselected -= 1
            if isempty(stops)
                delete!(unselected, template)
                push!(sweep.spare, stops)
            end
        end
    end
    # Activate blocks starting at this position
    while !isempty(pending) && (@inbounds pending[1])[1] <= position
        (_, stop, template) = heappop!(pending)
        heappush!(active, (stop, template))
        if template in sweep.selected
            sweep.selected_depth += 1
        else
            stops = get(unselected, template, nothing)
            if stops === nothing
                stops = isempty(sweep.spare) ? Int[] : pop!(sweep.spare)
                unselected[template] = stops
            end
            push!(stops, stop)
            sweep.n_unselected += 1
            heappush!(candidates, (stop, template), Reverse)
        end
    end
    # Select templates until the target depth is reached, preferring the template
    # whose block extends the furthest.
    # Every active block of an unselected template has an entry in `candidates`,
    # so it can't run empty while selected_depth < length(active).
    target = min(sweep.depth, length(active))
    while sweep.selected_depth < target
        (stop, template) = heappop!(candidates, Reverse)
        stop < position && continue
        stops = pop!(unselected, template, nothing)
        stops === nothing && continue
        push!(sweep.selected, template)
        sweep.selected_depth += length(stops)
        sweep.n_unselected -= length(stops)
        empty!(stops)
        push!(sweep.spare, stops)
    end
    # Remove stale entries, so `candidates` doesn't grow without bound in regions
    # where no selection is needed.
    if length(candidates) > 2 * sweep.n_unselected + 1024
        filter!(candidates) do (stop, template)
            stop >= position && haskey(unselected, template)
        end
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
