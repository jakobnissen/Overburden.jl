module Overburden

using XAM2: XAM2, BAMReader, BAMRecord, BAMWriter, SAMHeader, XAMError, cigar, get_header, is_mapped, left_position, query_name, reference_index, references
using BGZFLib: BGZFReader, BGZFError
using BinaryHeaps: heapify!, heappop!, heappush!
using Base.Order: By, Reverse
using CIGARStrings: CIGARStrings, BAMCIGAR, CIGARError, aligned_blocks
using OptParse: OptParse, @parser, arg, default, flag, integer, option, optparse, print_help, str

# Read at precompile time, so the executable doesn't need the Project.toml
const PROJECT_TOML = joinpath(dirname(@__DIR__), "Project.toml")
include_dependency(PROJECT_TOML)
const OVERBURDEN_VERSION = let
    m = match(r"^version\s*=\s*\"(.*)\"\s*$"m, read(PROJECT_TOML, String))
    m === nothing ? error("no version in Project.toml") : String(m[1])
end

"""
    OverburdenError(msg::String)

An error caused by bad user input, such as an unsorted BAM file.
The command-line interface reports it without a stack trace.
"""
struct OverburdenError <: Exception
    msg::String
end

Base.showerror(io::IO, e::OverburdenError) = print(io, e.msg)

"""
    SortChecker(n_references::Int)

Checks that records are sorted as the sweep needs: The records of each reference are
contiguous, and sorted by position. Unlike the SAM specification's coordinate order,
references may come in any order, e.g. from concatenating per-reference BAM files,
and records with no reference may be anywhere, since they don't take part in the sweep.
"""
mutable struct SortChecker
    # References whose records have been seen, i.e. which must not start again
    const seen::BitVector
    # Index of the reference of the last record with a reference, or 0 if none
    reference::UInt
    # Position of the last record with a reference, or 0 if it has no position
    position::UInt
end

SortChecker(n_references::Int) = SortChecker(falses(n_references), 0, 0)

"""
    check_sorted!(checker::SortChecker, record::BAMRecord, n::Int)::Nothing

Check that `record`, the `n`th record, is in order after the records previously
checked. Throw an `OverburdenError` if not.
"""
function check_sorted!(checker::SortChecker, record::BAMRecord, n::Int)::Nothing
    # Reference indices and positions are positive, so the conversions are lossless
    index = (@something reference_index(record) return nothing) % UInt
    position = something(left_position(record), 0) % UInt
    if index != checker.reference
        if @inbounds checker.seen[index]
            throw_unsorted(record, n, "the records of its reference are not contiguous")
        end
        @inbounds checker.seen[index] = true
        checker.reference = index
    elseif position < checker.position
        throw_unsorted(record, n, "its position precedes that of the record before it")
    end
    checker.position = position
    return nothing
end

@noinline function throw_unsorted(record::BAMRecord, n::Int, reason::String)
    name = something(query_name(record), "*")
    reference = something(XAM2.reference_name(record), "*")
    throw(OverburdenError("input is not coordinate sorted: record $n (\"$name\" on $reference): $reason"))
end

"""
    open_bam(path::String; n_workers::Int=Threads.nthreads(), copy::Bool=true)::BAMReader

Open the BAM file at `path`, and check that its header declares it coordinate sorted.
Throw an `OverburdenError` if the file is missing or not coordinate sorted, or an
`XAMError` or `BGZFError` if it is malformed.
"""
function open_bam(path::String; n_workers::Int = Threads.nthreads(), copy::Bool = true)
    isfile(path) || throw(OverburdenError("input file not found: " * path))
    io = open(path)
    bgzf = try
        BGZFReader(io; n_workers)
    catch
        close(io)
        rethrow()
    end
    reader = XAM2.try_init(BAMReader, bgzf; copy)
    reader isa XAMError && throw(reader)
    try
        check_sorted_header(reader)
    catch
        close(reader)
        rethrow()
    end
    return reader
end

function check_sorted_header(reader::BAMReader)::Nothing
    header = get_header(reader)
    if header isa XAMError
        throw(OverburdenError("malformed BAM header: " * sprint(showerror, header)))
    end
    hd = header.header
    if hd === nothing || hd.sort_order != XAM2.SortOrders.coordinate
        throw(OverburdenError("input is not coordinate sorted: header lacks \"@HD SO:coordinate\""))
    end
    return nothing
end

include("sweep.jl")
include("output.jl")
include("cli.jl")

end # module Overburden
