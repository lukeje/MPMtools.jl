module MRIutils

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
    ernst(α, TR, R₁)

Steady state signal from the Ernst equation for given flip angle α, TR, and R₁
    
α should be in radians, and time units of TR and R₁ should match.

    ernst(α, TR, T1=T₁)

Steady state signal using T₁ instead of R₁

α should be in radians, and time units of TR and T₁ should match.
"""
function ernst(α::Number, TR::Number, R₁::Number)
    sin(α) * (1 - exp(-TR * R₁)) / (1 - cos(α) * exp(-TR * R₁))
end

ernst(α, TR; T1::Number) = ernst(α, TR, one(T1)/T1)

"""
    ernstd(α, TR, R₁)

Steady state signal from the Ernst equation for given flip angle α, TR, and R₁.
    
α should be in degrees, and time units of TR and R₁ should match.

    ernstd(α, TR, T1=T₁)

Steady state signal using T₁ instead of R₁.

α should be in degrees, and time units of TR and R₁ should match.
"""
ernstd(α, TR, R₁) = ernst(deg2rad(α), TR, R₁)

ernstd(α, TR; T1::Number) = ernstd(α, TR, one(T1)/T1)


"""
    half_angle_tan(α)
Compute half angle tangent transform of α, where α is in radians.
Used to transform Ernst equation into analytically-soluble form.

# Reference
- Dathe and Helms, Phys. Med. Biol. (2010), "Exact algebraization of 
    the signal equation of spoiled gradient echo MRI".
    https://doi.org/10.1088/0031-9155/55/15/003
"""
half_angle_tan(α) = 2tan(0.5α)

"""
    optimalVFAangles(TR, R₁[, PDorR1=("PD" | "R1")])

Optimal flip angles for estimating PD or R₁ from dual flip angle VFA in radians.

# Reference
- Dathe and Helms, Phys. Med. Biol. (2010), "Exact algebraization of 
    the signal equation of spoiled gradient echo MRI".
    https://doi.org/10.1088/0031-9155/55/15/003 
"""
function optimalVFAangles(TR::Number, R₁::Number, PDorR1::String="R1")

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
    optimalVFAanglesd(TR, R₁[, PDorR1=("PD" | "R1")])

Optimal flip angles for estimating R₁ from VFA in degrees.

# Reference
- Dathe and Helms, Phys. Med. Biol. (2010), "Exact algebraization of 
    the signal equation of spoiled gradient echo MRI".
    https://doi.org/10.1088/0031-9155/55/15/003 
"""
function optimalVFAanglesd(TR, R₁, PDorR1::String="R1") 
    map(rad2deg, optimalVFAangles(TR, R₁, PDorR1))
end

"""
    dR1(SPD,ST1,dSPD,dST1,α_PD,α_T1,TRPD,TRT1)
Calculate propagation of uncertainty for R1 map.

In:
- SPD:  PDw signal at TE=0
- ST1:  T1w signal at TE=0
- dSPD: residual of mono-exponential fit of PDw signal
- dST1: residual of mono-exponential fit of T1w signal
- α_PD: flip angle of PDw signal
- α_T1: flip angle of T1w signal
- TRPD: repetition time of PDw signal
- TRT1: repetition time of T1w signal

Out: error for R1 in reciprocal units of TR units

# References
- https://en.wikipedia.org/wiki/Propagation_of_uncertainty
- Mohammadi et al. NeuroImage (2022), "Error quantification in
    multi-parameter mapping facilitates robust estimation and enhanced
    group level sensitivity"
    https://doi.org/10.1016/j.neuroimage.2022.119529
"""
function dR1(SPD,ST1,dSPD,dST1,α_PD,α_T1,TRPD,TRT1)
    
    τ_PD = half_angle_tan(α_PD)
    τ_T1 = half_angle_tan(α_T1)

    """
    Derivative of dual flip-angle R1 estimate with respect to first weighted 
    signal (S1). Because of symmetry in the R1 calculation, the derivative 
    with respect to the second weighted signal can be computed by permuting
    labels.
    """
    dR1_by_dS1(S1, S2, α1, α2, TR1, TR2) =
        (S1*α1/(2*TR1) - S2*α2/(2*TR2)) / (α1*(S1/α1 - S2/α2)^2) - 
        α1/(2*TR1*(S1/α1 - S2/α2))

    # dR1 calculation is symmetric with respect to the two weighted contrasts
    sqrt( dR1_by_dS1(SPD,ST1,τ_PD,τ_T1,TRPD,TRT1)^2*dSPD^2 + 
          dR1_by_dS1(ST1,SPD,τ_T1,τ_PD,TRT1,TRPD)^2*dST1^2 )
    
end

"""
    dPD(SPD,ST1,dSPD,dST1,α_PD,α_T1,TRPD,TRT1)
Calculate propagation of uncertainty for PD map.

In:
- SPD:  PDw signal at TE=0
- ST1:  T1w signal at TE=0
- dSPD: residual of mono-exponential fit of PDw signal
- dST1: residual of mono-exponential fit of T1w signal
- α_PD: flip angle of PDw signal
- α_T1: flip angle of T1w signal
- TRPD: repetition time of PDw signal
- TRT1: repetition time of T1w signal

Out: error for A in arbitrary units (a.u.)
 
# References
- https://en.wikipedia.org/wiki/Propagation_of_uncertainty
- Mohammadi et al. NeuroImage (2022), "Error quantification in 
    multi-parameter mapping facilitates robust estimation and enhanced 
    group level sensitivity." 
    https://doi.org/10.1016/j.neuroimage.2022.119529
"""
function dPD(SPD,ST1,dSPD,dST1,α_PD,α_T1,TRPD,TRT1)

    τ_PD = half_angle_tan(α_PD)
    τ_T1 = half_angle_tan(α_T1)

    """
    Derivative of dual flip-angle A (PD) estimate with respect to first 
    weighted signal (S1). Because of symmetry in the R1 calculation, the 
    derivative with respect to the second weighted signal can be computed by 
    permuting labels.
    """
    dPD_by_dS1(S1,S2,α1,α2,TR1,TR2) =
        S1*S2*TR2*α1*(TR1*α2/α1 - TR2*α1/α2) / (S1*TR2*α1 - S2*TR1*α2)^2 -
            S2*(TR1*α2/α1 - TR2*α1/α2)/(S1*TR2*α1 - S2*TR1*α2)

    # dPD calculation is symmetric with respect to the two weighted contrasts
    sqrt( dPD_by_dS1(SPD,ST1,τ_PD,τ_T1,TRPD,TRT1)^2*dSPD^2 +
          dPD_by_dS1(ST1,SPD,τ_T1,τ_PD,TRT1,TRPD)^2*dST1^2 )

end

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
