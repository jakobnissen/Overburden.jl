# Tests of pass 2 and the entry point. Included from runtests.jl.

using Overburden: overburden
using XAM2: BAMReader, flags, left_position, query_name, reference_name
using MemoryViews: MemoryView

# (name, flag, reference, position) of each record of the BAM file at `path`
function read_summary(path::String)
    return BAMReader(open(path); n_workers = 1) do reader
        map(reader) do record
            (
                String(something(query_name(record), "*")), flags(record),
                String(something(reference_name(record), "*")),
                something(left_position(record), 0),
            )
        end
    end
end

# Text of the header of the BAM file at `path`
function header_text(path::String)::String
    return BAMReader(open(path); n_workers = 1) do reader
        String(MemoryView(get_header(reader)))
    end
end

@testset "Output" begin
    rng = Xoshiro(4)
    for _ in 1:10
        records = random_records(rng, rand(rng, 1:300))
        input = write_records(TEST_DIR, records)
        for depth in (1, 3, 20)
            output = tempname(TEST_DIR) * ".bam"
            selected = select_templates(input, depth; n_workers = 1)
            (; n_templates, n_written) = overburden(input, output, depth; n_workers = 1)
            @test n_templates == length(selected)
            # The records of the selected templates, in input order
            expected = filter(r -> template_hash(r[1]) in selected, read_summary(input))
            @test read_summary(output) == expected
            @test n_written == length(expected)
        end
    end
end

@testset "Output header" begin
    header = HEADER * "@PG\tID:overburden\tPN:overburden\n"
    input = sam_to_bam(TEST_DIR, header * SORTED)
    output = tempname(TEST_DIR) * ".bam"
    overburden(input, output, 1; n_workers = 1, command_line = "overburden\tx\ny")
    lines = split(chomp(header_text(output)), '\n')
    @test lines[1:4] == split(chomp(header), '\n')
    @test lines[5] == "@PG\tID:overburden.1\tPN:overburden\tCL:overburden x y\tPP:overburden\tVN:$(Overburden.OVERBURDEN_VERSION)"
    @test length(lines) == 5
end

@testset "Output errors" begin
    input = sam_to_bam(TEST_DIR, HEADER * SORTED)
    @test occursin("is the input", error_message(() -> overburden(input, input, 1)))
    # Also through a symlink
    link = tempname(TEST_DIR)
    symlink(input, link)
    @test occursin("is the input", error_message(() -> overburden(input, link, 1)))
    missing_dir = joinpath(TEST_DIR, "missing", "out.bam")
    @test occursin("directory not found", error_message(() -> overburden(input, missing_dir, 1)))
    # An existing output is not replaced, even if it can't be written
    for mode in (0o644, 0o444)
        existing = tempname(TEST_DIR)
        write(existing, "old")
        chmod(existing, mode)
        @test occursin("already exists", error_message(() -> overburden(input, existing, 1; n_workers = 1)))
        @test read(existing, String) == "old"
    end
    # A dangling symlink is an existing output too: creating it would create its target
    dangling = tempname(TEST_DIR)
    target = joinpath(TEST_DIR, "nowhere.bam")
    symlink(target, dangling)
    @test occursin("already exists", error_message(() -> overburden(input, dangling, 1; n_workers = 1)))
    @test islink(dangling) && !ispath(target)
end

@testset "Output created early" begin
    input = sam_to_bam(TEST_DIR, HEADER * SORTED)
    # The output is created before the first pass, so this is a fresh file during the run
    output = tempname(TEST_DIR) * ".bam"
    overburden(input, output, 1; n_workers = 1)
    @test length(read_summary(output)) == 3
    # A failing first pass removes the output, so it doesn't block the next run
    unsorted = sam_to_bam(TEST_DIR, HEADER * sam_line("a", 0, "chr1", 50, "10M") * sam_line("b", 0, "chr1", 10, "10M"))
    unsorted_output = tempname(TEST_DIR) * ".bam"
    @test occursin("not coordinate sorted", error_message(() -> overburden(unsorted, unsorted_output, 1; n_workers = 1)))
    @test !ispath(unsorted_output)
    # The output can't be created in a directory that isn't writable, and nothing is replaced
    dir = joinpath(TEST_DIR, "readonly_dir")
    mkdir(dir)
    chmod(dir, 0o555)
    try
        @test_throws Base.IOError overburden(input, joinpath(dir, "out.bam"), 1; n_workers = 1)
        @test isempty(readdir(dir))
    finally
        chmod(dir, 0o755)
    end
end

@testset "Entry point" begin
    input = sam_to_bam(TEST_DIR, HEADER * SORTED)
    output = tempname(TEST_DIR) * ".bam"
    fresh_output() = tempname(TEST_DIR) * ".bam"
    # Run the entry point, and return its exit code, stdout and stderr
    function run_main(args::Vector{String})
        out = tempname(TEST_DIR)
        err = tempname(TEST_DIR)
        code = redirect_stdio(() -> Overburden.main(args); stdout = out, stderr = err)
        return (code, read(out, String), read(err, String))
    end
    (code, out, err) = run_main([input, output])
    @test code == 0
    @test isempty(out)
    @test err == "Wrote 4 records of 4 templates\n"
    @test length(read_summary(output)) == 4
    @test run_main(["-d", "1", input, fresh_output()])[1] == 0
    @test run_main([input, fresh_output(), "--depth=1"])[1] == 0
    # The output must not exist yet
    (code, out, err) = run_main([input, output])
    @test code == 1
    @test occursin("already exists", err)
    for args in [["-h"], ["--help"], [input, "--help"]]
        (code, out, err) = run_main(args)
        @test code == 0
        @test occursin("Usage: overburden", out)
        @test isempty(err)
    end
    for args in [["-V"], ["--version"]]
        @test run_main(args) == (0, "overburden $(Overburden.OVERBURDEN_VERSION)\n", "")
    end
    # Invalid arguments
    (code, out, err) = run_main(String[])
    @test code == 2
    @test occursin("Usage: overburden", err)
    for (args, msg) in [
            ([input], "too few arguments"),
            ([input, output, "x"], ""),
            (["--bogus", input, output], "--bogus"),
            (["-d", "0", input, output], "minimum"),
            (["-d", "x", input, output], "integer"),
        ]
        (code, out, err) = run_main(args)
        @test code == 2
        @test isempty(out)
        @test occursin(msg, err)
    end
    # Errors
    for (args, msg) in [
            ([joinpath(TEST_DIR, "missing.bam"), output], "not found"),
            ([input, input], "is the input"),
        ]
        (code, out, err) = run_main(args)
        @test code == 1
        @test occursin("ERROR: ", err) && occursin(msg, err)
    end
end
