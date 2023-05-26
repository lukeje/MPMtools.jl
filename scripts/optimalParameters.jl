using ArgParse
using Optim

include("./MPMtools.jl/MRIutils.jl")
import .MRIutils as MPM

function parse_commandline()
    s = ArgParseSettings("Calculate optimal flip angles (α1 and α1) and repetition times (TR1 and TR2) for dual flip angle R1 and PD measurements.")

    @add_arg_table! s begin
        "nlines"
            help = "number of acquired lines"
            required = true
            arg_type = Float64
        "scantime"
            help = "acquisition time to fill (s)"
            required = true
            arg_type = Float64
        "R1"
            help = "representative R1 rate (1/s)"
            required = true
            arg_type = Float64
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
    end

    add_arg_group!(s, "whether to optimise for \"PD\", \"R1\", or \"both\"", exclusive=true, required=true)
    @add_arg_table! s begin
        "--onlyPD"
            action = :store_true
        "--onlyR1"
            action = :store_true
        "--both"
            action = :store_true
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()

    R1     = parsed_args["R1"]
    TRmin  = parsed_args["TRmin"]
    FAmax  = deg2rad(parsed_args["FAmax"])
    TRsum  = parsed_args["scantime"]/parsed_args["nlines"]
    
    @assert (TRsum > TRmin) "The requested scantime and number of lines is not consistent with the minimal TR. Please relax your input parameters and try again."

    constrainedTR2(TR1, s) = (TRsum - TR1 - TRmin)*s + TRmin # ensure TR2 ∈ [TRmin, TRsum - TRmin] given fit parameters s ∈ [0,1] and TR1 ∈ [TRmin, TRsum - TRmin]

    # dPD and dR1 have same argument list, so define them here rather than repeating them below
    dargs(x) = ( MPM.ernst(x[1], x[2], R1), MPM.ernst(x[3], constrainedTR2(x[2], x[4]), R1), 1.0, 1.0, x[1], x[3], x[2], constrainedTR2(x[2], x[4]) )

    # choose the fitting function based on which quantitative parameter the output parameters should be optimal for estimating
    if parsed_args["onlyPD"]
        fitfun = x -> MPM.dPD(dargs(x)...)
        initialoptimum = "PD"
    elseif parsed_args["onlyR1"]
        fitfun = x -> MPM.dR1(dargs(x)...)
        initialoptimum = "R1"
    elseif parsed_args["both"]
        fitfun = x -> MPM.dPD(dargs(x)...) + MPM.dR1(dargs(x)...)/R1
        initialoptimum = "PD"
    end

    # start from optimum for equal TRs
    angles = MPM.optimalVFAangles(TRsum/2, R1, initialoptimum)
    angles = min.(angles, FAmax-0.01) # ensure angles within bounds

    opt = optimize(fitfun, [0.0, TRmin, 0.0, 0.0], [FAmax, TRsum-TRmin, FAmax, 1.0], [angles[1], TRsum/2, angles[2], 1.0 - 0.01])
    xopt = opt.minimizer

    # precision in output
    rounddigit(x) = round(x; digits=2)

    println("α1:  $(rounddigit(rad2deg(xopt[1])))°")
    println("TR1: $(rounddigit(1e3*xopt[2])) ms")
    println("α2:  $(rounddigit(rad2deg(xopt[3])))°")
    println("TR2: $(rounddigit(1e3*(constrainedTR2(xopt[2], xopt[4])))) ms")

end

main()
