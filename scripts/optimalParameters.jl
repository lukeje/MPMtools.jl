using ArgParse
using Optim

include("./MPMtools.jl/MRIutils.jl")
import .MRIutils as MPM

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--TRmin"
            help = "minimum TR for readout and excitation (s)"
            default = Inf64
            arg_type = Float64
        "--nlines"
            help = "number of acquired lines"
            required = true
            arg_type = Float64
        "--R1"
            help = "representative R1 rate (1/s)"
            required = true
            arg_type = Float64
        "--scantime"
            help = "acquisition time to fill (s)"
            required = true
            arg_type = Float64
        "--PDorR1"
            help = "whether to optimise for \"PD\", \"R1\", or \"both\""
            required = true
            arg_type = String
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()

    R1 = parsed_args["R1"]
    TRmin = parsed_args["TRmin"]
    TRsum = parsed_args["scantime"]/(2*parsed_args["nlines"])

    PDorR1 = parsed_args["PDorR1"]
    if PDorR1 == "PD"
        fitfun = x -> MPM.dPD(MPM.ernst(x[1], x[2], R1), MPM.ernst(x[3], TRsum-x[2], R1), 1.0, 1.0, x[1], x[3], x[2], TRsum-x[2])
        initialoptimum = PDorR1
    elseif PDorR1 == "R1"
        fitfun = x -> MPM.dR1(MPM.ernst(x[1], x[2], R1), MPM.ernst(x[3], TRsum-x[2], R1), 1.0, 1.0, x[1], x[3], x[2], TRsum-x[2])
        initialoptimum = PDorR1
    elseif PDorR1 == "both"
        fitfun = x -> MPM.dPD(MPM.ernst(x[1], x[2], R1), MPM.ernst(x[3], TRsum-x[2], R1), 1.0, 1.0, x[1], x[3], x[2], TRsum-x[2]) +
                    MPM.dR1(MPM.ernst(x[1], x[2], R1), MPM.ernst(x[3], TRsum-x[2], R1), 1.0, 1.0, x[1], x[3], x[2], TRsum-x[2])/R1
        initialoptimum = "PD"
    else
        error("could not process $PDorR1")
    end

    angles = MPM.optimalVFAangles(TRsum/2, R1, initialoptimum)

    opt = optimize(fitfun, [0, 0, 0], [3π/2, TRsum, 3π/2], [angles[1], TRsum/2, angles[2]])
    xopt = opt.minimizer

    print("α1: $(rad2deg(xopt[1]))\n")
    print("TR1: $(1e3*xopt[2])\n")
    print("α2: $(rad2deg(xopt[3]))\n")
    print("TR2: $(1e3*(TRsum - xopt[2]))\n")

    @assert all([xopt[2], (TRsum - xopt[2])] .≥ TRmin) "One or both TRs were less than TRmin. Please relax your input parameters and try again."

end

main()