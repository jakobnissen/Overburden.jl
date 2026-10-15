# Command-line interface

# Help is written by hand, since OptParse's `print_help` can't be trimmed (OptParse 0.4.2)
const HELP = """
Write a near-minimal subset of the templates of a coordinate-sorted BAM file, keeping
the depth at every position at least D, or the input's depth where that is lower.
Templates are written in full or not at all.

Usage: overburden [OPTIONS] <INPUT> <OUTPUT>

Arguments:
  <INPUT>          Coordinate-sorted input BAM file
  <OUTPUT>         Output BAM file, which must not exist

Options:
  -d, --depth <D>  Minimum depth to keep [default: 20]
  -h, --help       Print this help and exit
  -V, --version    Print the version and exit
"""

@parser struct Options
    depth = default(option("-d", "--depth", integer("D"; min = 1)), 20)
    # Handled before parsing, but included so they're not rejected
    help = flag("-h", "--help")
    version = flag("-V", "--version")
    input = arg(str("INPUT"))
    output = arg(str("OUTPUT"))
end

const ERROR_TYPES = Union{OverburdenError, XAMError, BGZFError, SystemError, Base.IOError}

function print_error(e::ERROR_TYPES)
    print(Core.stderr, "ERROR: ")
    # Some of Base's `showerror` methods can't be trimmed
    if e isa SystemError
        print(Core.stderr, e.prefix, ": ", Libc.strerror(e.errnum))
    elseif e isa Base.IOError
        print(Core.stderr, e.msg)
    else
        showerror(Core.stderr, e)
    end
    println(Core.stderr)
    return nothing
end

"""
    main(args::Vector{String})::Cint

Entry point of the `overburden` executable. Return the exit code: 0 on success,
1 on errors, and 2 on invalid command-line arguments.
"""
function (@main)(args::Vector{String})::Cint
    if isempty(args)
        print(Core.stderr, HELP)
        return 2
    end
    # Help and version are handled before parsing, since the positional arguments
    # are not required with them.
    for arg in args
        arg == "--" && break
        if arg in ("-h", "--help")
            print(Core.stdout, HELP)
            return 0
        elseif arg in ("-V", "--version")
            println(Core.stdout, "overburden ", OVERBURDEN_VERSION)
            return 0
        end
    end
    options = optparse(Options, args)
    options === nothing && return 2
    command_line = join(pushfirst!(copy(args), "overburden"), ' ')
    return try
        (; n_templates, n_written) = overburden(
            options.input, options.output, options.depth; command_line,
        )
        println(Core.stderr, "Wrote ", n_written, " records of ", n_templates, " templates")
        0
    catch e
        e isa ERROR_TYPES || rethrow()
        print_error(e)
        1
    end
end
