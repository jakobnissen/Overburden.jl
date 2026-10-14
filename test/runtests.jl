using Test
using Overburden: Overburden, OverburdenError, select_templates
using XAM2: XAM2, BAMRecord, BAMWriter, SAMReader, XAMError, get_header

const HEADER = """
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:1000
@SQ\tSN:chr2\tLN:500
"""

# Write the SAM text `sam` as a BAM file in `dir`, and return its path
function sam_to_bam(dir::String, sam::String)
    path = tempname(dir) * ".bam"
    SAMReader(IOBuffer(sam)) do reader
        BAMWriter(open(path, "w"), get_header(reader); n_workers = 1) do writer
            for record in reader
                bam = BAMRecord(writer)
                copy!(bam, record)
                write(writer, bam)
            end
        end
    end
    return path
end

sam_line(name, flag, rname, pos, cigar) =
    "$name\t$flag\t$rname\t$pos\t60\t$cigar\t*\t0\t0\t*\t*\n"

const SORTED = join(
    [
        sam_line("a", 0, "chr1", 5, "4M"),
        sam_line("b", 0, "chr1", 5, "10M"),
        sam_line("c", 0, "chr1", 900, "3M2N3M"),
        sam_line("d", 0, "chr2", 1, "4M"),
        # Unmapped, but placed at its mate's position
        sam_line("e", 4, "chr2", 20, "*"),
        sam_line("f", 4, "*", 0, "*"),
        sam_line("g", 4, "*", 0, "*"),
    ]
)

# Return the message of the `OverburdenError` thrown by `f()`
function error_message(f)
    try
        f()
    catch e
        e isa OverburdenError && return e.msg
        rethrow()
    end
    error("no error thrown")
end

# Number of templates selected from the BAM file at `path`
count(path; depth = 20) = length(select_templates(path, depth; n_workers = 1))

mktempdir() do dir
    @testset "Sorted input" begin
        # e, f and g are unmapped
        @test count(sam_to_bam(dir, HEADER * SORTED)) == 4
        @test count(sam_to_bam(dir, HEADER)) == 0
    end

    @testset "Header must be coordinate sorted" begin
        for header in [
                "@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:500\n",
                replace(HEADER, "coordinate" => "unsorted"),
                replace(HEADER, "coordinate" => "queryname"),
                replace(HEADER, "\tSO:coordinate" => ""),
            ]
            path = sam_to_bam(dir, header * SORTED)
            @test occursin("SO:coordinate", error_message(() -> count(path)))
        end
    end

    @testset "Records must be sorted within each reference" begin
        for (records, reason) in [
                (
                    [sam_line("a", 0, "chr1", 10, "4M"), sam_line("b", 0, "chr1", 9, "4M")],
                    "position precedes",
                ),
                # chr1 starts again after chr2
                (
                    [
                        sam_line("a", 0, "chr1", 1, "4M"),
                        sam_line("b", 0, "chr2", 1, "4M"),
                        sam_line("c", 0, "chr1", 9, "4M"),
                    ],
                    "not contiguous",
                ),
                # Records with no reference don't reset the position
                (
                    [
                        sam_line("a", 0, "chr1", 10, "4M"),
                        sam_line("b", 4, "*", 0, "*"),
                        sam_line("c", 0, "chr1", 9, "4M"),
                    ],
                    "position precedes",
                ),
            ]
            path = sam_to_bam(dir, HEADER * join(records))
            msg = error_message(() -> count(path))
            @test occursin("not coordinate sorted", msg)
            @test occursin("(\"$(records[end][1])\" on chr1)", msg)
            @test occursin(reason, msg)
        end
    end

    @testset "References may be in any order" begin
        records = [
            # Records with no reference may be anywhere
            sam_line("a", 4, "*", 0, "*"),
            sam_line("b", 0, "chr2", 5, "4M"),
            sam_line("c", 0, "chr2", 7, "4M"),
            sam_line("d", 4, "*", 0, "*"),
            sam_line("e", 0, "chr1", 1, "4M"),
            sam_line("f", 0, "chr1", 9, "4M"),
        ]
        @test count(sam_to_bam(dir, HEADER * join(records))) == 4
    end

    @testset "Bad files" begin
        @test occursin("not found", error_message(() -> count(joinpath(dir, "missing.bam"))))
        # A SAM file is not BGZF compressed
        sam_path = joinpath(dir, "input.sam")
        write(sam_path, HEADER * SORTED)
        @test_throws Exception count(sam_path)
    end

    @testset "Entry point" begin
        path = sam_to_bam(dir, HEADER * SORTED)
        # Silence the output of the entry point
        redirect(f) = redirect_stdio(f; stdout = devnull, stderr = devnull)
        @test redirect(() -> Overburden.main([path])) == 0
        @test redirect(() -> Overburden.main(String[])) == 2
        @test redirect(() -> Overburden.main([joinpath(dir, "missing.bam")])) == 1
    end
end

include("sweep.jl")
