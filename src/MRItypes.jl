module MRItypes

export WeightedContrast, WeightedMultiechoContrast

"""
    WeightedContrast(signal, flipangle, TR[, TE])

Composite type to store variable flip angle (VFA) MRI data
"""
struct WeightedContrast
    signal
    flipangle
    TR
    TE::Union{<:Number,Missing}
    τ

    # define explicit inner constructor as τ is computed from the input
    # and we need to check some parameter bounds
    function WeightedContrast(signal, flipangle, TR, TE)

        @assert (TR > 0) "TR must be greater than zero!"
        @assert (TE ≥ 0) | ismissing(TE) "TE must be greater than or equal to zero or missing!"

        # populate τ field (half angle tangent transform) using provided flipangle
        τ = half_angle_tan(flipangle)

        new(signal, flipangle, TR, TE, τ)
    end
end

# outer constructor: TE is missing if not provided
WeightedContrast(signal, flipangle, TR) = WeightedContrast(signal, flipangle, TR, missing)


"""
    WeightedMultiechoContrast([w(TE₁), w(TE₂), ...])

Vector type to store multiecho variable flip angle (VFA) MRI data
"""
struct WeightedMultiechoContrast
    signal::Vector{<:Number}
    flipangle
    TR
    TE::Vector{<:Number}
    τ

    # use constructor to check consistency of data
    function WeightedMultiechoContrast(echoList::Vector{WeightedContrast})
        @assert length(echoList) > 1 "WeightedMultiechoContrast requires more than one element!"
        
        signal = [w.signal for w in echoList]
        flipangle = echoList[1].flipangle
        TR = echoList[1].TR
        TE = [w.TE for w in echoList]
        τ = echoList[1].τ

        for w in echoList
            @assert w.flipangle == flipangle "All flip angles must match!"
            @assert w.TR == TR "All TRs must match!"
            @assert w.TE > 0 "All TEs must be greater than zero!"
        end

        @assert (length(unique(TE)) > 1) "There must be more than one unique TE!"

        new(signal,flipangle,TR,TE,τ)
    end
    
end


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

end