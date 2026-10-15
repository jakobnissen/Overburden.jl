# Pass 2: Write the records of the selected templates.

"""
    write_templates(writer::BAMWriter, reader::BAMReader, selected::Set{UInt64})::Int

Read all remaining records of `reader`, and write those whose [`template_hash`](@ref) is
in `selected` to `writer`, copying their raw bytes. Return the number of records written.
Throw an `XAMError` if the file is malformed.
"""
function write_templates(writer::BAMWriter, reader::BAMReader, selected::Set{UInt64})::Int
    n_written = 0
    while true
        record = XAM2.try_read!(reader)
        record === nothing && break
        record isa XAMError && throw(record)
        if template_hash(query_name(record)) in selected
            write(writer, record)
            n_written += 1
        end
    end
    return n_written
end

"""
    output_header(reader::BAMReader, command_line::Union{String, Nothing})::SAMHeader

The header of `reader`, with an added `@PG` line for Overburden.
"""
function output_header(reader::BAMReader, command_line::Union{String, Nothing})::SAMHeader
    # `open_bam` has checked the header is not an `XAMError`
    header = get_header(reader)::SAMHeader
    # Header fields can't contain tabs, newlines etc.
    cl = if command_line === nothing
        nothing
    else
        map(c -> c in ('\t', '\n', '\r', '\0') ? ' ' : c, command_line)
    end
    return XAM2.add_program(header, "overburden"; pn = "overburden", vn = OVERBURDEN_VERSION, cl)
end

"""
    create_output(input::String, output::String)::Nothing

Check that `output` can be created, and create it as an empty file. Throw an
`OverburdenError` if `output` is the file `input`, if it already exists (including as a
symbolic link, which `touch` would follow), or if its directory does not exist. Throw a
`Base.IOError` if creating the file fails. The check and the creation are not atomic, so a
file created by another process in between is not detected.
"""
function create_output(input::String, output::String)::Nothing
    if ispath(output) && samefile(input, output)
        throw(OverburdenError("output file is the input file: " * output))
    end
    (ispath(output) || islink(output)) && throw(OverburdenError("output file already exists: " * output))
    output_dir = dirname(abspath(output))
    isdir(output_dir) || throw(OverburdenError("output directory not found: " * output_dir))
    touch(output)
    return nothing
end

# Level 5 rather than BGZFLib's 6, since compression is the bottleneck of pass 2, and
# level 5 makes it about 15% faster, for about 2% more output.
const DEFAULT_COMPRESS_LEVEL = 5

"""
    overburden(
        input::String, output::String, depth::Int;
        n_workers::Int = Threads.nthreads(), compress_level::Int = $(DEFAULT_COMPRESS_LEVEL),
        command_line::Union{String, Nothing} = nothing,
    )::@NamedTuple{n_templates::Int, n_written::Int}

Read the coordinate-sorted BAM file `input`, and write the records of a near-minimal set
of its templates to the BAM file `output`, such that at every reference position `p`,
the output's depth is at least `min(depth, d(p))`, where `d(p)` is the input's depth.
See [`select_templates`](@ref).

The output is compressed with BGZF at `compress_level`, from 1 (fastest) to 12
(smallest). It has the input's header with an added `@PG` line, whose `CL` field is
`command_line`, if given. Return the number of templates selected, and the number of
records written.

Throw an `OverburdenError` on bad input, such as unsorted records, or if `output` already
exists, is the same file as `input`, or its directory does not exist. The output is created
before the first pass, so an existing file is never replaced. If an error is thrown after
it is created, it is removed.
"""
function overburden(
        input::String,
        output::String,
        depth::Int;
        n_workers::Int = Threads.nthreads(),
        compress_level::Int = DEFAULT_COMPRESS_LEVEL,
        command_line::Union{String, Nothing} = nothing,
    )::@NamedTuple{n_templates::Int, n_written::Int}
    in(compress_level, 1:12) || throw(ArgumentError("compress_level must be in 1:12"))
    reader = open_bam(input; n_workers, copy = false)
    return try
        create_output(input, output)
        try
            selected = select_templates(reader, depth)
            header = output_header(reader, command_line)
            XAM2.seek_first_record(reader)
            n_written = BAMWriter(open(output, "w"), header, references(reader); n_workers, compress_level) do writer
                write_templates(writer, reader, selected)
            end
            (; n_templates = length(selected), n_written)
        catch
            # The output was created above, so it's ours to remove. Don't leave an empty or
            # truncated file, which would block the next run.
            rm(output; force = true)
            rethrow()
        end
    finally
        close(reader)
    end
end
