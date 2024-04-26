using Pkg
Pkg.activate(dirname(@__DIR__))
using ArgParse
using MPMtools.MRIutils: optimalVFAparameters

function parse_commandline()
    s = ArgParseSettings("Calculate optimal flip angles (α1 and α1) and repetition times (TR1 and TR2) for dual flip angle R1 and PD measurements.")

    s.epilog = """
               examples: \n
               1. Given a matrix size of 100×100, a maximum scan time of 12 minutes, and a typical R1 of 1.25 / s,
               the optimal TR and flip angle pair can be computed using:\n
               \ua0\ua0julia $(s.prog) \$((100*100)) \$((12*60)) 1.25 --onlyR1\n
               \n\n
               2. If we know the minimum and maximum TR values rather than having a maximum scan time, 
               then we can specify a minimum TR (here 18.7 ms), 
               the sum of the TR values as the scan time (here 18.7 ms + 23.7 ms),
               and set the number of lines to 1. 
               Here we also specify a maximum flip angle of 20°:\n
               \ua0\ua0julia $(s.prog) 1 0.0424 1.25 \\\n
               \ua0\ua0\ua0\ua0--TRmin 0.0187 --onlyR1 --FAmax 20\n
               This case reproduces the result from Weiskopf, et al. (2013, https://doi.org/10.3389/fnins.2013.00095).
               """

    @add_arg_table! s begin
        "scantime"
            help = "acquisition time to fill (s)"
            required = true
            arg_type = Float64
        "R1"
            help = "representative R1 rate (1/s)"
            required = true
            arg_type = Float64
            nargs    = '*'
        "--nlines"
            help = "number of acquired lines"
            required = false
            arg_type = Int64
            default  = 1
        "--nvolumes"
            help = "number of volumes to use"
            required = false
            arg_type = Int64
            default  = 2
        "--TRmin"
            help = "minimum TR for readout and excitation (s)"
            required = false
            arg_type = Float64
            default  = 0.0
        "--FAmax"
            help = "maximum flip angle for excitation (°)"
            required = false
            arg_type = Float64
            default  = 270.0
        "--B1range"
            help = "minimum, interval, and maximum of points in B1 range to optimise over in fractions of nominal flip angle"
            required = false
            arg_type = Float64
            nargs    = 3
            default  = [1.0,1.0,1.0]
    end

    add_arg_group!(s, "whether to optimise for \"PD\", \"R1\", or \"both\". \"both\" requires an argument in [0,1] specifying the relative weighting of \"PD\" and \"R1\"", exclusive=true, required=true)
    @add_arg_table! s begin
        "--onlyPD"
            action = :store_true
        "--onlyR1"
            action = :store_true
        "--both"
            arg_type = Float64
            default  = 0.5
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()

    R1     = parsed_args["R1"]
    nvols  = parsed_args["nvolumes"]
    TRmin  = parsed_args["TRmin"]
    FAmax  = deg2rad(parsed_args["FAmax"])
    TRsum  = parsed_args["scantime"]/parsed_args["nlines"]
    B1     = range(start=parsed_args["B1range"][1],step=parsed_args["B1range"][2],stop=parsed_args["B1range"][3])

    if parsed_args["onlyPD"]
        PDorR1 = "PD"
    elseif parsed_args["onlyR1"]
        PDorR1 = "R1"
    else 
        PDorR1 = parsed_args["both"]
    end
    
    @assert (TRsum > 2TRmin) "The requested scantime and number of lines is not consistent with the minimal TR. Please relax your input parameters and try again."

    FA,TR = optimalVFAparameters(TRsum, R1, nvols, B1, PDorR1=PDorR1, TRmin=TRmin, FAmax=FAmax)

    # precision in output
    rounddigit(x) = round(x; digits=2)

    for n in 1:nvols
        println("α$(n):  $(rounddigit(rad2deg(FA[n])))°")
        println("TR$(n): $(rounddigit(1e3*TR[n])) ms")
    end

end

main()
