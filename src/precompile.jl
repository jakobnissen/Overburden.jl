# Precompile the command-line interface to reduce its latency.
# Run the entry point quietly, so the workload stays representative if it's refactored.

@setup_workload begin
    data = joinpath(dirname(@__DIR__), "test", "data")
    @compile_workload begin
        mktempdir() do dir
            run_cli(["-d", "1", joinpath(data, "sorted.bam"), joinpath(dir, "out.bam")], true)
            # Error path: Unsorted input
            run_cli([joinpath(data, "unsorted.bam"), joinpath(dir, "unsorted.bam")], true)
        end
        run_cli(["--help"], true)
        precompile(main, (Vector{String},))
    end
end
