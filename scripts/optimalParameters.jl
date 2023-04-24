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
    TRsum  = parsed_args["scantime"]/(2*parsed_args["nlines"])

    if parsed_args["onlyPD"]
        fitfun = x -> MPM.dPD(MPM.ernst(x[1], x[2], R1), MPM.ernst(x[3], TRsum-x[2], R1), 1.0, 1.0, x[1], x[3], x[2], TRsum-x[2])
        initialoptimum = "PD"
    elseif parsed_args["onlyR1"]
        fitfun = x -> MPM.dR1(MPM.ernst(x[1], x[2], R1), MPM.ernst(x[3], TRsum-x[2], R1), 1.0, 1.0, x[1], x[3], x[2], TRsum-x[2])
        initialoptimum = "R1"
    elseif parsed_args["both"]
        fitfun = x -> MPM.dPD(MPM.ernst(x[1], x[2], R1), MPM.ernst(x[3], TRsum-x[2], R1), 1.0, 1.0, x[1], x[3], x[2], TRsum-x[2]) +
                      MPM.dR1(MPM.ernst(x[1], x[2], R1), MPM.ernst(x[3], TRsum-x[2], R1), 1.0, 1.0, x[1], x[3], x[2], TRsum-x[2])/R1
        initialoptimum = "PD"
    end

    angles = MPM.optimalVFAangles(TRsum/2, R1, initialoptimum)

    opt = optimize(fitfun, [0, 0, 0], [3π/2, TRsum, 3π/2], [angles[1], TRsum/2, angles[2]])
    xopt = opt.minimizer

    print("α1:  $(rad2deg(xopt[1]))°\n")
    print("TR1: $(1e3*xopt[2]) ms\n")
    print("α2:  $(rad2deg(xopt[3]))°\n")
    print("TR2: $(1e3*(TRsum - xopt[2])) ms\n")

    @assert all([xopt[2], (TRsum - xopt[2])] .≥ TRmin) "One or both TRs were less than TRmin. Please relax your input parameters and try again."

end

main()