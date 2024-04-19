module MRIutils

using Optim
using LinearAlgebra
using ..MRItypes
using ..MRItypes: half_angle_tan

"""
    ernstangle(TR, R₁)

Angle maximising the Ernst equation in radians for given TR and R₁.

# Reference
- Equation (11.52), BKZ, "Handbook of MRI Pulse Sequences" (2004)
"""
ernstangle(TR::Number, R₁::Number) = acos(exp(-TR * R₁))

"""
    ernstangled(TR, R₁)

Angle maximising the Ernst equation in degrees for given TR and R₁.

# Reference
- Equation (11.52), BKZ, "Handbook of MRI Pulse Sequences" (2004)
"""
ernstangled(TR, R₁) = rad2deg(ernstangle(TR, R₁))


"""
    ernst(α, TR, R₁[, PD=PD])

Steady state signal from the Ernst equation for given flip angle α, TR, and R₁.
Optionally scale the signal by PD.
    
α should be in radians, and time units of TR and R₁ should match.

    ernst(α, TR, T1=T₁[, PD=PD])

Steady state signal using T₁ instead of R₁

α should be in radians, and time units of TR and T₁ should match.
"""
function ernst(α::Number, TR::Number, R₁::Number; PD::Number=one(TR))
    signal = PD * sin(α) * (1 - exp(-TR * R₁)) / (1 - cos(α) * exp(-TR * R₁))
    return WeightedContrast(signal, α, TR)
end

ernst(α, TR; T1::Number, PD=one(TR)) = ernst(α, TR, inv(T1), PD=PD)

"""
    ernstd(α, TR, R₁)

Steady state signal from the Ernst equation for given flip angle α, TR, and R₁.
    
α should be in degrees, and time units of TR and R₁ should match.

    ernstd(α, TR, T1=T₁)

Steady state signal using T₁ instead of R₁.

α should be in degrees, and time units of TR and R₁ should match.
"""
ernstd(α, TR, R₁; PD=one(TR)) = ernst(deg2rad(α), TR, R₁, PD=PD)

ernstd(α, TR; T1::Number, PD=one(TR)) = ernstd(α, TR, one(T1)/T1, PD=PD)

"""
    exponentialDecay(S::WeightedContrast, lambda, TElist)

Scale signal of WeightedContrast by exponential decay S.signal = exp(-lambda * TE) for each TE in TElist.
Returns a WeightedMultiechoContrast type.
"""
function exponentialDecay(S::WeightedContrast, lambda, TElist)
    WeightedMultiechoContrast([WeightedContrast(S.signal*exp(-lambda * te), S.flipangle , S.TR, te) for te in TElist])
end


"""
    dR1(SPD,ST1,dSPD,dST1)
Calculate propagation of uncertainty for R1 map.

In:
- SPD:  PDw WeightedContrast at TE=0
- ST1:  T1w WeightedContrast at TE=0
- dSPD: residual of mono-exponential fit of PDw signal
- dST1: residual of mono-exponential fit of T1w signal

Out: error for R1 in reciprocal units of TR units

    dR1(R₁,PD,dSPD,dST1)
Calculate propagation of uncertainty for R1 map using synthetic signal values for given R₁.

In:
- α_PD: flip angle of PDw signal in radians
- α_T1: flip angle of T1w signal in radians
- TRPD: repetition time of PDw signal
- TRT1: repetition time of T1w signal
- dSPD: variance of PDw signal
- dST1: variance of T1w signal

Out: error for R1 in reciprocal units of TR units

# References
- https://en.wikipedia.org/wiki/Propagation_of_uncertainty
- Mohammadi et al. NeuroImage (2022), "Error quantification in
    multi-parameter mapping facilitates robust estimation and enhanced
    group level sensitivity"
    https://doi.org/10.1016/j.neuroimage.2022.119529
"""
function dR1(SPD::WeightedContrast,ST1::WeightedContrast,dSPD::Number,dST1::Number)

    """
    Derivative of dual flip-angle R1 estimate with respect to first weighted 
    signal (S1). Because of symmetry in the R1 calculation, the derivative 
    with respect to the second weighted signal can be computed by permuting
    labels.
    """
    dR1_by_dS1(S1, S2, α1, α2, TR1, TR2) =
        (S1*α1/(2*TR1) - S2*α2/(2*TR2)) / (α1*(S1/α1 - S2/α2)^2) - α1/(2*TR1*(S1/α1 - S2/α2))

    # dR1 calculation is symmetric with respect to the two weighted contrasts
    sqrt( dR1_by_dS1(SPD.signal,ST1.signal,SPD.τ,ST1.τ,SPD.TR,ST1.TR)^2*dSPD^2 + 
          dR1_by_dS1(ST1.signal,SPD.signal,ST1.τ,SPD.τ,ST1.TR,SPD.TR)^2*dST1^2 )
    
end

dR1(R₁,dSPD,dST1,α_PD,α_T1,TRPD,TRT1) = dR1(ernst(α_PD,TRPD,R₁),ernst(α_T1,TRT1,R₁),dSPD,dST1)


"""
    dPD(SPD,ST1,dSPD,dST1)
Calculate propagation of uncertainty for PD map.

In:
- SPD:  PDw signal at TE=0
- ST1:  T1w signal at TE=0
- dSPD: residual of mono-exponential fit of PDw signal
- dST1: residual of mono-exponential fit of T1w signal

Out: error for A in arbitrary units (a.u.)

    dPD(R₁,dSPD,dST1)
Calculate propagation of uncertainty for PD map using synthetic signal values for given R₁.

In:
- α_PD: flip angle of PDw signal in radians
- α_T1: flip angle of T1w signal in radians
- TRPD: repetition time of PDw signal
- TRT1: repetition time of T1w signal
- dSPD: variance of PDw signal
- dST1: variance of T1w signal

Out: error for A in arbitrary units (a.u.)
 
# References
- https://en.wikipedia.org/wiki/Propagation_of_uncertainty
- Mohammadi et al. NeuroImage (2022), "Error quantification in 
    multi-parameter mapping facilitates robust estimation and enhanced 
    group level sensitivity." 
    https://doi.org/10.1016/j.neuroimage.2022.119529
"""
function dPD(SPD::WeightedContrast,ST1::WeightedContrast,dSPD::Number,dST1::Number)

    """
    Derivative of dual flip-angle A (PD) estimate with respect to first 
    weighted signal (S1). Because of symmetry in the R1 calculation, the 
    derivative with respect to the second weighted signal can be computed by 
    permuting labels.
    """
    dPD_by_dS1(S1,S2,α1,α2,TR1,TR2) =
        S1*S2*TR2*α1*(TR1*α2/α1 - TR2*α1/α2) / (S1*TR2*α1 - S2*TR1*α2)^2 - S2*(TR1*α2/α1 - TR2*α1/α2)/(S1*TR2*α1 - S2*TR1*α2)

    # dPD calculation is symmetric with respect to the two weighted contrasts
    sqrt( dPD_by_dS1(SPD.signal,ST1.signal,SPD.τ,ST1.τ,SPD.TR,ST1.TR)^2*dSPD^2 + 
          dPD_by_dS1(ST1.signal,SPD.signal,ST1.τ,SPD.τ,ST1.TR,SPD.TR)^2*dST1^2 )

end

dPD(R₁,dSPD,dST1,α_PD,α_T1,TRPD,TRT1) = dPD(ernst(α_PD,TRPD,R₁),ernst(α_T1,TRT1,R₁),dSPD,dST1)


"""
    dT1([SPD,ST1,...], [dSPD,dST1,...])
Calculate propagation of uncertainty for T1 and A map.

In:
- SPD:  PDw WeightedContrast at TE=0
- ST1:  T1w WeightedContrast at TE=0
- dSPD: residual of mono-exponential fit of PDw signal
- dST1: residual of mono-exponential fit of T1w signal

Out: error for A and R1 in reciprocal units of TR units

    dT1(R₁,[dSPD,dST1,...],[αPD,αT1,...],,[TRPD,TRT1,...])
Calculate propagation of uncertainty for T1 and A map using synthetic signal values for given R₁.

In:
- α_PD: flip angle of PDw signal in radians
- α_T1: flip angle of T1w signal in radians
- TRPD: repetition time of PDw signal
- TRT1: repetition time of T1w signal
- dSPD: variance of PDw signal
- dST1: variance of T1w signal

Out: error for A in arbitrary units and T1 in TR units

# References
- https://en.wikipedia.org/wiki/Propagation_of_uncertainty
- Mohammadi et al. NeuroImage (2022), "Error quantification in
    multi-parameter mapping facilitates robust estimation and enhanced
    group level sensitivity"
    https://doi.org/10.1016/j.neuroimage.2022.119529
- Helms et al. Magn. Reson. Med. (2011), "Identification of signal bias 
    in the variable flip angle method by linear display of the 
    algebraic ernst equation", 
    [doi:10.1002/mrm.22849](https://doi.org/10.1002/mrm.22849)
    
"""
function dT1(W::Vector{WeightedContrast},dS::Vector{<:Number})

    S  = [w.signal for w in W]
    τ  = [w.τ      for w in W]
    TR = [w.TR     for w in W]

    y = S./τ
    D = hcat(ones(length(S)), -S.*τ./(2TR))

    dA  = zero(eltype(S))
    dT1 = zero(eltype(S))
    for n in 1:length(S)
        y′ = zero(y)
        y′[n] = one(S[n])/τ[n]

        D′ = zero(D)
        D′[n,2] = -one(S[n])*τ[n]/(2TR[n])

        σ = D\(y′ .- D′*(D\y))

        dA  += first(σ)^2 * dS[n]^2
        dT1 += last(σ)^2  * dS[n]^2
    end
    
    return (sqrt(dA),sqrt(dT1)) 
end

dT1(R₁,dS,α,TR) = dT1([ernst(α,TR,R₁) for (α,TR) in zip(α,TR)],dS)


"""
    optimalDFAangles(TR, R₁[, PDorR1=("R1" | "PD")])

Optimal flip angles for estimating R₁ or PD from dual flip angle R1 mapping in radians.

# Reference
- Dathe and Helms, Phys. Med. Biol. (2010), "Exact algebraization of 
    the signal equation of spoiled gradient echo MRI".
    https://doi.org/10.1088/0031-9155/55/15/003 
"""
function optimalDFAangles(TR::Number, R₁::Number, PDorR1::String="R1")

    if PDorR1 == "PD"
        # Equation (27)
        scaling = #(0.49030, 3.14611)
            (0.4903044753219954, 
             3.1461052210182947)
    elseif PDorR1 == "R1"
        # Equation (21)
        scaling = (sqrt(2)-1, sqrt(2)+1)
    else
        error("PDorR1 must be either \"PD\" or \"R1\"")
    end

    # Equations (21) and (27)
    map(s -> 2atan(s * 0.5half_angle_tan(ernstangle(TR, R₁))), scaling)
end

"""
    optimalDFAanglesd(TR, R₁[, PDorR1=("R1" | "PD")])

Optimal flip angles for estimating R₁ or PD from dual flip angle R1 mapping in degrees.

# Reference
- Dathe and Helms, Phys. Med. Biol. (2010), "Exact algebraization of 
    the signal equation of spoiled gradient echo MRI".
    https://doi.org/10.1088/0031-9155/55/15/003 
"""
function optimalDFAanglesd(TR, R₁, PDorR1::String="R1") 
    map(rad2deg, optimalDFAangles(TR, R₁, PDorR1))
end


"""
    α1, α2, TR1, TR2 = optimalDFAparameters(TRsum, R₁[, PDorR1=("R1" | "PD" | x ∈ [0,1]), TRmin=0.0, FAmax=3π/2])

Optimal repetition times and flip angles for estimating R₁ and/or PD from dual flip angle R1 mapping.

# Notes
- Units of TRsum, R₁, and TRmin must be consistent. Output TRs will be in the same units.
- All angles are in radians.

# Reference
- TBC
"""
function optimalDFAparameters(TRsum, R₁; PDorR1::Union{String,Number}="R1", TRmin=0.0, FAmax=3π/2)
    
    @assert (TRsum > 2TRmin) "The requested TRsum is not consistent with the minimal TR. Please relax your input parameters and try again."

    # computes allowed TR2 ∈ [TRmin, TRsum - TRmin] given fit parameters s ∈ [0,1] and TR1 ∈ [TRmin, TRsum - TRmin]
    constrainedTR2(TR1, s) = (TRsum - TR1 - TRmin)*s + TRmin

    # dPD and dR1 have same argument list, so define them here rather than repeating them below
    dargs(x) = ( R₁, 1.0, 1.0, x[1], x[2], x[3], constrainedTR2(x[3], x[4]) )

    # choose the fitting function based on which quantitative parameter the output parameters should be optimal for estimating
    if PDorR1=="PD"
        PDR1fraction = 0.0
    elseif PDorR1=="R1"
        PDR1fraction = 1.0
    elseif isa(PDorR1,Number) && (0.0 <= PDorR1 <= 1.0)
        PDR1fraction = PDorR1
    else
        error("PDorR1 must be either \"PD\", \"R1\", or a relative weighting in [0,1]. Was $(PDorR1).")
    end

    if PDR1fraction==0.0
        fitfun = x -> dPD(dargs(x)...)^2
        initialoptimum = "PD"
    elseif PDR1fraction==1.0
        fitfun = x -> dR1(dargs(x)...)^2
        initialoptimum = "R1"
    else
        fitfun = x -> (1.0 - PDR1fraction)*dPD(dargs(x)...)^2 + PDR1fraction*(dR1(dargs(x)...)/R₁)^2
        initialoptimum = "PD"
    end

    # start from optimum for equal TRs
    # -0.01 so that optimiser does not start on edge of allowed range
    angles = optimalDFAangles(TRsum/2, R₁, initialoptimum)
    angles = min.(angles, FAmax - 0.01) # enforce FAmax in initial conditions
    x0(angles) = [angles..., TRsum/2, 1.0 - 0.01]

    # check both orders of initial flip angle guesses
    opt = [optimize(fitfun, [0.0, 0.0, TRmin, 0.0], [FAmax, FAmax, TRsum-TRmin, 1.0], x0(a)) for a in (angles, reverse(angles))]
    xopt = opt[argmin(o.minimum for o in opt)].minimizer

    # α1, α2, TR1, TR2
    return xopt[1], xopt[2], xopt[3], constrainedTR2(xopt[3], xopt[4])
end

function optimalDFAparameters(TR1, TR2, R₁; PDorR1::Union{String,Number}="R1", FAmax=3π/2)

    # dPD and dR1 have same argument list, so define them here rather than repeating them below
    dargs(x) = ( R₁, 1.0, 1.0, x[1], x[2], TR1, TR2 )

    # choose the fitting function based on which quantitative parameter the output parameters should be optimal for estimating
    if PDorR1=="PD"
        PDR1fraction = 0.0
    elseif PDorR1=="R1"
        PDR1fraction = 1.0
    elseif isa(PDorR1,Number) && (0.0 <= PDorR1 <= 1.0)
        PDR1fraction = PDorR1
    else
        error("PDorR1 must be either \"PD\", \"R1\", or a relative weighting in [0,1]. Was $(PDorR1).")
    end

    if PDR1fraction==0.0
        fitfun = x -> dPD(dargs(x)...)^2
        initialoptimum = "PD"
    elseif PDR1fraction==1.0
        fitfun = x -> dR1(dargs(x)...)^2
        initialoptimum = "R1"
    else
        fitfun = x -> (1.0 - PDR1fraction)*dPD(dargs(x)...)^2 + PDR1fraction*(dR1(dargs(x)...)/R₁)^2
        initialoptimum = "PD"
    end

    # start from optimum for equal TRs
    # -0.01 so that optimiser does not start on edge of allowed range
    angles = optimalDFAangles((TR1+TR2)/2, R₁, initialoptimum)
    angles = min.(angles, FAmax - 0.01) # enforce FAmax in initial conditions
    x0 = [angles[1], angles[2]]

    # check both orders of initial flip angle guesses
    opt = [optimize(fitfun, [0.0, 0.0], [FAmax, FAmax], x) for x in (x0, reverse(x0))]
    xopt = opt[argmin(o.minimum for o in opt)].minimizer

    # α1, α2, TR1, TR2
    return xopt[1], xopt[2], TR1, TR2
end

"""
    α, TR = optimalVFAparameters(TRsum, R₁, nvolumes[, PDorR1=("R1" | "PD" | x ∈ [0,1]), TRmin=0.0, FAmax=3π/2])

Optimal repetition times and flip angles for estimating R₁ and/or PD from dual flip angle R1 mapping.

# Notes
- Units of TRsum, R₁, and TRmin must be consistent. Output TRs will be in the same units.
- All angles are in radians.

# Reference
- TBC
"""
function optimalVFAparameters(TRsum, R₁, nvolumes, B1; PDorR1::Union{String,Number}="R1", TRmin=0.0, FAmax=3π/2)
    
    @assert (TRsum > nvolumes*TRmin) "The requested TRsum is not consistent with the minimal TR. Please relax your input parameters and try again."

    # computes allowed TR2 ∈ [TRmin, TRsum - TRmin] given fit parameters s ∈ [0,1] and TR1 ∈ [TRmin, TRsum - TRmin]
    function constrainedTR2(TR1, S)
        TR = [TR1]
        sizehint!(TR,length(S)+1)
        tsum = TR1
        for (i,s) in pairs(S)
            push!(TR,(TRsum - tsum - (length(S)-i+1)*TRmin)*s + TRmin)
            tsum += last(TR)
        end
        return TR
    end

    # choose the fitting function based on which quantitative parameter the output parameters should be optimal for estimating
    if PDorR1=="PD"
        PDT1fraction = 0.0
    elseif PDorR1=="R1"
        PDT1fraction = 1.0
    elseif isa(PDorR1,Number) && (0.0 <= PDorR1 <= 1.0)
        PDT1fraction = PDorR1
    else
        error("PDorR1 must be either \"PD\", \"R1\", or a relative weighting in [0,1]. Was $(PDorR1).")
    end

    function dT1_2(α,TR)
        d = zeros(Float64,2)
        for b in B1 
            d .+= dT1(R₁, ones(Float64,nvolumes), b.*α, TR).^2
        end
        return d
    end

    if PDT1fraction==0.0
        fitfun = x -> first(dT1_2(x[begin:nvolumes], constrainedTR2(x[nvolumes+1],x[nvolumes+2:end])))
        initialoptimum = "PD"
    elseif PDT1fraction==1.0
        fitfun = x -> last( dT1_2(x[begin:nvolumes], constrainedTR2(x[nvolumes+1],x[nvolumes+2:end])))
        initialoptimum = "R1"
    else
        fitfun = x -> dot([one(PDT1fraction) - PDT1fraction, PDT1fraction*(R₁^2)], 
                            dT1_2(x[begin:nvolumes], constrainedTR2(x[nvolumes+1],x[nvolumes+2:end])))
        initialoptimum = "PD"
    end

    # start from optimum for equal TRs
    # -0.01 so that optimiser does not start on edge of allowed range
    angles = range(optimalDFAangles(TRsum/nvolumes, R₁, initialoptimum)...,nvolumes)
    angles = min.(angles, FAmax - 0.01) # enforce FAmax in initial conditions
    x0(angles) = vcat(angles, TRsum/nvolumes, ones(Float64,nvolumes-1) .- 0.01)

    lower = vcat(zeros(Float64,nvolumes),      TRmin,                    zeros(Float64,nvolumes-1))
    upper = vcat(ones(Float64,nvolumes)*FAmax, TRsum-(nvolumes-1)*TRmin, ones(Float64,nvolumes-1))

    opt = optimize(fitfun, lower, upper, x0(angles))
    xopt = opt.minimizer
    
    # [α1, α2, ...], [TR1, TR2, ...]
    return xopt[begin:nvolumes], constrainedTR2(xopt[nvolumes+1],xopt[nvolumes+2:end])
end

optimalVFAparameters(TRsum, R₁, nvolumes; PDorR1::Union{String,Number}="R1", TRmin=0.0, FAmax=3π/2) = 
    optimalVFAparameters(TRsum, R₁, nvolumes, one(R₁):one(R₁); PDorR1=PDorR1, TRmin=TRmin, FAmax=FAmax)

"""
    inversionRecovery(R₁, TI, η)

The long-TR inversion recovery signal for given R₁, inversion time TI and inversion efficiency η

    inversionRecovery(R₁, TI)

The long-TR inversion recovery signal for given R₁ and inversion time TI assuming perfect inversion (η = 1)
"""
function inversionRecovery(R₁::Number, TI::Number, η::Number)

    @assert (0 < η ≤ 1) "Inversion efficiency η must be between 0 and 1"
    
    return 1 - 2η * exp(-TI * R₁)
end

inversionRecovery(R₁, TI) = inversionRecovery(R₁, TI, 1)


end
