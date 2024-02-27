using Pkg
Pkg.activate(dirname(@__DIR__))
using ArgParse
using MPMtools.MRIutils: optimalDFAparameters

function parse_commandline()
    s = ArgParseSettings("Calculate optimal flip angles (α1 and α1) for dual flip angle R1 and PD measurements given repetition times and R1.")

    s.epilog = """
               examples: \n
               ...
               """

    @add_arg_table! s begin
        "TR1"
            help = "first TR (s)"
            required = true
            arg_type = Float64
        "TR2"
            help = "second TR (s)"
            required = true
            arg_type = Float64
        "R1"
            help = "representative R1 rate (1/s)"
            required = true
            arg_type = Float64
        "--FAmax"
            help = "maximum flip angle for excitation (°)"
            required = false
            arg_type = Float64
            default  = 270.0
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
    FAmax  = deg2rad(parsed_args["FAmax"])
    TR1    = parsed_args["TR1"]
    TR2    = parsed_args["TR2"]

    if parsed_args["onlyPD"]
        PDorR1 = "PD"
    elseif parsed_args["onlyR1"]
        PDorR1 = "R1"
    else 
        PDorR1 = parsed_args["both"]
    end
    
    xopt = optimalDFAparameters(TR1, TR2, R1, PDorR1=PDorR1, FAmax=FAmax)

    # precision in output
    rounddigit(x) = round(x; digits=2)

    println("α1:  $(rounddigit(rad2deg(xopt[1])))°")
    println("TR1: $(rounddigit(1e3*xopt[3])) ms")
    println("α2:  $(rounddigit(rad2deg(xopt[2])))°")
    println("TR2: $(rounddigit(1e3*xopt[4])) ms")

end

main()
