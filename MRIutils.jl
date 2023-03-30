module MRIutils

"""
    ernstangle(TR, R₁)

Angle maximising the Ernst equation in radians for given TR and R₁
"""
function ernstangle(TR::Number, R₁::Number)
    # Equation (11.52), BKZ, "Handbook of MRI Pulse Sequences" (2004)
    acos(exp(-TR * R₁))
end

"""
    ernstangled(TR, R₁)

Angle maximising the Ernst equation in degrees for given TR and R₁
"""
ernstangled(TR, R₁) = rad2deg(ernstangle(TR, R₁))

"""
    ernst(α, TR, R₁)

Steady state signal from the Ernst equation for given flip angle α, TR, and R₁
    
α should be in radians, and time units of TR and R₁ should match.
"""
function ernst(α::Number, TR::Number, R₁::Number)
    sin(α) * (1 - exp(-TR * R₁)) / (1 - cos(α) * exp(-TR * R₁))
end

"""
    ernstd(α, TR, R₁)

Steady state signal from the Ernst equation for given flip angle α, TR, and R₁
    
α should be in degrees, and time units of TR and R₁ should match.
"""
ernstd(α, TR, R₁) = ernst(deg2rad(α), TR, R₁)

"""
    optimalVFAangles(TR, R₁[, PDorR1=("PD" | "R1")])

Optimal flip angles for estimating PD or R₁ from dual flip angle VFA in radians.

Uses the solution from Dathe and Helms, "Exact algebraization of the signal equation of spoiled gradient echo MRI", 
Physics in Medicine and Biology (2010) [doi:10.1088/0031-9155/55/15/003](https://dx.doi.org/10.1088/0031-9155/55/15/003) 
"""
function optimalVFAangles(TR::Number, R₁::Number, PDorR1::String="R1")

    if PDorR1 ==  "PD"
        # Equation (27)
        scaling = (0.49030, 3.14611)
    elseif PDorR1 == "R1"
        # Equation (21)
        scaling = (sqrt(2)-1, sqrt(2)+1)
    else
        error("PDorR1 must be either \"PD\" or \"R1\"")
    end

    # Equations (21) and (27)
    return map(s -> 2atan(s * tan(0.5*ernstangle(TR, R₁))), scaling)
end

"""
    optimalVFAanglesd(TR, R₁[, PDorR1=("PD" | "R1")])

Optimal flip angles for estimating R₁ from VFA in degrees

Uses the solution from Dathe and Helms, "Exact algebraization of the signal equation of spoiled gradient echo MRI", 
Physics in Medicine and Biology (2010) [doi:10.1088/0031-9155/55/15/003](https://dx.doi.org/10.1088/0031-9155/55/15/003) 
"""
function optimalVFAanglesd(TR, R₁, PDorR1::String="R1") 
    map(t -> rad2deg(t), optimalVFAangles(TR, R₁, PDorR1))
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
