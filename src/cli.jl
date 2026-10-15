# Command-line interface

@parser struct Options
    @description """
    Write a near-minimal subset of the templates of a coordinate-sorted BAM file, keeping
    the depth at every position at least D, or the input's depth where that is lower.
    Templates are written in full or not at all."""
    "Minimum depth to keep"
    depth = default(option("-d", "--depth", integer("D"; min = 1)), 20)
    # Handled before parsing, but included so they're not rejected, and appear in the help
    "Print this help and exit"
    help = flag("-h", "--help")
    "Print the version and exit"
    version = flag("-V", "--version")
    "Coordinate-sorted input BAM file"
    input = arg(str("INPUT"))
    "Output BAM file, which must not exist"
    output = arg(str("OUTPUT"))
end

print_options_help(io::IO) = print_help(io, OptParse.lift(Options), String[]; progname = "overburden")

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
        print_options_help(Core.stderr)
        return 2
    end
    # Help and version are handled before parsing, since the positional arguments
    # are not required with them.
    for arg in args
        arg == "--" && break
        if arg in ("-h", "--help")
            print_options_help(Core.stdout)
            return 0
        elseif arg in ("-V", "--version")
            println(Core.stdout, "overburden ", OVERBURDEN_VERSION)
            return 0
        end
    end
    options = optparse(Options, args; progname = "overburden")
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
