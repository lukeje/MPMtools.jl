module MRIutils

"""
    ernstangle(TR, R₁)

The solution of the Ernst equation in radians for given TR and R₁
"""
function ernstangle(TR, R₁)
    # Equation (11.52), BKZ, "Handbook of MRI Pulse Sequences" (2004)
    acos(exp(-TR * R₁))
end

"""
    ernstangled(TR, R₁)

The solution of the Ernst equation in degrees for given TR and R₁
"""
function ernstangled(TR, R₁)
    rad2deg(ernstangle(TR, R₁))
end

"""
    ernst(α, TR, R₁)

The steady state signals from the Ernst equation for given flip angle α, TR, and R₁
    
α should be in radians, and time units of TR and R₁ should match.
"""
function ernst(α, TR, R₁)
    sin(α) * (1 - exp(-TR * R₁)) / (1 - cos(α) * exp(-TR * R₁))
end

"""
    ernstd(α, TR, R₁)

The ssteady state signals from the Ernst equation for given flip angle α, TR, and R₁
    
α should be in degrees, and time units of TR and R₁ should match.
"""
function ernstd(α, TR, R₁)
    ernst(deg2rad(α), TR, R₁)
end

"""
    optimalVFAangles(TR, R₁[, PDorR1=("PD" | "R1")])

Optimal flip angles for estimating R₁ from VFA in radians.

Uses the solution from Dathe and Helms, "Exact algebraization of the signal equation of spoiled gradient echo MRI", 
Physics in Medicine and Biology (2010) [doi:10.1088/0031-9155/55/15/003](https://dx.doi.org/10.1088/0031-9155/55/15/003) 
"""
function optimalVFAangles(TR, R₁, PDorR1::String="R1")

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
    return 2 * atan(scaling * tan(ernstangle(TR, R₁) / 2))
end

"""
    optimalVFAanglesd(TR, R₁[, PDorR1=("PD" | "R1")])

Optimal flip angles for estimating R₁ from VFA in degrees

Uses the solution from Dathe and Helms, "Exact algebraization of the signal equation of spoiled gradient echo MRI", 
Physics in Medicine and Biology (2010) [doi:10.1088/0031-9155/55/15/003](https://dx.doi.org/10.1088/0031-9155/55/15/003) 
"""
function optimalVFAanglesd(TR, R₁, PDorR1::String="R1")
    rad2deg(optimalVFAangles(TR, R₁, PDorR1))
end

"""
    inversionRecovery(R₁, TI, η)

The long-TR inversion recovery signal for given R₁, inversion time TI and inversion efficiency η
"""
function inversionRecovery(R₁, TI, η)

    @assert (0 < η ≤ 1) "Inversion efficiency η must be between 0 and 1"
    
    return 1 - 2η * exp(-TI * R₁)
end

end