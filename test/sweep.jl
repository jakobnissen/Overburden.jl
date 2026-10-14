# Tests of the pass-1 selection. Included from runtests.jl.

using Random: Random, Xoshiro
using Overburden: template_hash

# A record of a synthetic test file
struct TestRecord
    name::String
    flag::Int
    reference::String # "*" if none
    position::Int # 0 if none
    cigar::String # "*" if none
end

# Directory for the test BAM files, deleted at exit
const TEST_DIR = mktempdir()

const REFERENCE_LENGTHS = Dict("chr1" => 1000, "chr2" => 500)

is_unmapped(r::TestRecord) = (r.flag & 4) != 0

# The positions covered by `record`, computed independently of `aligned_blocks`
function covered_positions(record::TestRecord)::Vector{Int}
    covered = Int[]
    (is_unmapped(record) || record.cigar == "*") && return covered
    p = record.position
    for m in eachmatch(r"(\d+)([MIDNSHP=X])", record.cigar)
        len = parse(Int, m[1])
        op = only(m[2])
        if op in "MD=X"
            append!(covered, p:(p + len - 1))
        end
        if op in "MD=XN"
            p += len
        end
    end
    return covered
end

# Depth per (reference, position) of the records for which `keep(record)` is true
function depths(keep, records::Vector{TestRecord})::Dict{Tuple{String, Int}, Int}
    result = Dict{Tuple{String, Int}, Int}()
    for record in records
        keep(record) || continue
        for p in covered_positions(record)
            key = (record.reference, p)
            result[key] = get(result, key, 0) + 1
        end
    end
    return result
end

# Sort as in a coordinate sorted file, write as BAM, and return the path
function write_records(dir::String, records::Vector{TestRecord})::String
    order = Dict("chr1" => 1, "chr2" => 2, "*" => 3)
    sorted = sort(records; by = r -> (order[r.reference], r.position))
    sam = join(sam_line(r.name, r.flag, r.reference, r.position, r.cigar) for r in sorted)
    return sam_to_bam(dir, HEADER * sam)
end

# Check that the selection from `records` with the given depth satisfies the
# requirements, and return the number of selected templates.
function check_selection(dir::String, records::Vector{TestRecord}, depth::Int)::Int
    path = write_records(dir, records)
    selected = select_templates(path, depth; n_workers = 1)
    is_selected(r::TestRecord) = template_hash(r.name) in selected
    full = depths(Returns(true), records)
    out = depths(is_selected, records)
    # Never below the target depth
    @test all(get(out, key, 0) >= min(depth, d) for (key, d) in full)
    # Templates covering nothing are never selected.
    covering = Set(template_hash(r.name) for r in records if !isempty(covered_positions(r)))
    @test selected ⊆ covering
    # No singletons: All records of a template are selected, or none.
    names = unique(r.name for r in records)
    for name in names
        segments = filter(r -> r.name == name, records)
        @test allequal(is_selected, segments)
    end
    return length(selected)
end

function random_cigar(rng)::String
    parts = String[]
    rand(rng) < 0.2 && push!(parts, "$(rand(rng, 1:5))S")
    push!(parts, "$(rand(rng, 1:40))M")
    for _ in 1:rand(rng, 0:2)
        op = rand(rng, ("I", "D", "N", "=", "X"))
        len = op == "N" ? rand(rng, 1:60) : rand(rng, 1:5)
        push!(parts, "$len$op", "$(rand(rng, 1:30))M")
    end
    rand(rng, Bool) && push!(parts, "$(rand(rng, 1:5))S")
    return join(parts)
end

# Random templates of one to three segments. Most segments are mapped, some are
# unmapped (placed at a mapped mate, or with no reference), and a few are mapped with
# no CIGAR. Flags other than unmapped don't matter, so are not set.
function random_records(rng, n_templates::Int)::Vector{TestRecord}
    records = TestRecord[]
    for i in 1:n_templates
        name = "t$i"
        mate = nothing
        for _ in 1:rand(rng, 1:3)
            if rand(rng) < 0.85
                reference = rand(rng, ("chr1", "chr2"))
                position = rand(rng, 1:(REFERENCE_LENGTHS[reference] - 150))
                cigar = rand(rng) < 0.02 ? "*" : random_cigar(rng)
                record = TestRecord(name, 0, reference, position, cigar)
                mate = record
            elseif mate !== nothing && rand(rng, Bool)
                record = TestRecord(name, 4, mate.reference, mate.position, "*")
            else
                record = TestRecord(name, 4, "*", 0, "*")
            end
            push!(records, record)
        end
    end
    return records
end

@testset "Template hash" begin
    @test template_hash(nothing) == template_hash("*")
    names = ["", "a", "b", "ab", "ba", "abcdefgh", "abcdefghi", "abcdefgh\0"]
    append!(names, ["read$i" for i in 1:1000])
    @test allunique(template_hash.(names))
end

@testset "Depth targets" begin
    # Positions at or below the depth keep all templates
    records = [
        TestRecord("a", 0, "chr1", 1, "10M"),
        TestRecord("b", 0, "chr1", 5, "10M"),
        TestRecord("c", 0, "chr1", 30, "5M"),
    ]
    @test check_selection(TEST_DIR, records, 2) == 3
    @test check_selection(TEST_DIR, records, 1) == 3
    # One long read covers the whole region
    records = [
        TestRecord("long", 0, "chr1", 1, "100M"),
        [TestRecord("s$i", 0, "chr1", i, "10M") for i in 1:10:91]...,
    ]
    @test check_selection(TEST_DIR, records, 1) == 1
    # Overlapping mates count twice
    records = [
        TestRecord("pair", 0, "chr1", 1, "10M"),
        TestRecord("pair", 0, "chr1", 1, "10M"),
    ]
    @test check_selection(TEST_DIR, records, 2) == 1
    # Splice gaps don't cover
    records = [
        TestRecord("spliced", 0, "chr1", 1, "10M20N10M"),
        TestRecord("gap", 0, "chr1", 15, "5M"),
    ]
    @test check_selection(TEST_DIR, records, 1) == 2
    # Selection is per template, also across references.
    records = [
        TestRecord("x", 0, "chr1", 1, "10M"),
        TestRecord("x", 0, "chr2", 1, "10M"),
        TestRecord("y", 0, "chr2", 1, "10M"),
        TestRecord("x", 4, "*", 0, "*"),
    ]
    @test check_selection(TEST_DIR, records, 1) == 1
    # The sweep resets between references, even when an unmapped record with a
    # reference comes between them.
    records = [
        TestRecord("a", 0, "chr1", 100, "10M"),
        TestRecord("b", 4, "chr2", 1, "*"),
        TestRecord("c", 0, "chr2", 105, "10M"),
    ]
    @test check_selection(TEST_DIR, records, 1) == 2
end

@testset "Random templates" begin
    rng = Xoshiro(1)
    for _ in 1:20
        records = random_records(rng, rand(rng, 1:300))
        n_covering = length(Set(r.name for r in records if !isempty(covered_positions(r))))
        for depth in (1, 2, 5, 20)
            n = check_selection(TEST_DIR, records, depth)
            @test n <= n_covering
        end
        # With a depth above the maximum, all templates must be selected
        @test check_selection(TEST_DIR, records, 10_000) == n_covering
    end
end

# The minimum number of reads to select from `records` to satisfy the depth target,
# by brute force over all subsets.
function minimum_selection(records::Vector{TestRecord}, depth::Int)::Int
    covered = map(r -> Set(covered_positions(r)), records)
    full = depths(Returns(true), records)
    best = length(records)
    for mask in 0:(2^length(records) - 1)
        n = count_ones(mask)
        n < best || continue
        ok = all(full) do ((_, p), d)
            out = sum(i -> (mask >> (i - 1)) & 1 == 1 && p in covered[i], eachindex(records))
            out >= min(depth, d)
        end
        ok && (best = n)
    end
    return best
end

@testset "Optimal on contiguous single-segment reads" begin
    rng = Xoshiro(2)
    for _ in 1:40
        records = [
            TestRecord("r$i", 0, "chr1", rand(rng, 1:30), "$(rand(rng, 1:12))M")
                for i in 1:rand(rng, 1:11)
        ]
        depth = rand(rng, 1:3)
        @test check_selection(TEST_DIR, records, depth) == minimum_selection(records, depth)
    end
end

@testset "High depth" begin
    # Enough stale candidates to trigger compaction of the candidate heap
    rng = Xoshiro(3)
    records = random_records(rng, 4000)
    n_covering = length(Set(r.name for r in records if !isempty(covered_positions(r))))
    for depth in (1, 5, 30)
        @test check_selection(TEST_DIR, records, depth) < n_covering
    end
end
