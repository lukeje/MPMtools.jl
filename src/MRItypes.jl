module MRItypes

export WeightedContrast, WeightedMultiechoContrast, SESTEcontrast

"""
    WeightedContrast(signal, flipangle, TR[, TE, TM=0])

Composite type to store variable flip angle (VFA) MRI data
"""
struct WeightedContrast
    signal
    flipangle
    TR
    TE::Union{<:Number,Missing}
    τ
    TM

    # define explicit inner constructor as τ is computed from the input
    # and we need to check some parameter bounds
    function WeightedContrast(signal, flipangle, TR, TE, TM=zero(TR))

        @assert                  TR > zero(TR) "TR must be greater than zero!"
        @assert ismissing(TE) || TE ≥ zero(TE) "TE must be greater than or equal to zero or missing!"
        @assert                  TM ≥ zero(TR) "TM must be greater than zero!"

        # populate τ field (half angle tangent transform) using provided flipangle
        τ = half_angle_tan(flipangle)

        new(signal, flipangle, TR, TE, τ, TM)
    end
end

# outer constructor: TE missing if not provided
WeightedContrast(signal, flipangle, TR) = WeightedContrast(signal, flipangle, TR, missing)

"""
    SESTEcontrast(signal, flipangle, TR[, TE, TM=0])

Composite type to store spin echo (SE)/stimulated echo (STE) MRI data
"""
struct SESTEcontrast
    SE::WeightedContrast
    STE::WeightedContrast

    function SESTEcontrast(SE,STE)
        # check contrasts are consistent
        # we do not check TE, as this might be different or not depending on definition
        @assert size(SE.signal) == size(STE.signal)
        @assert SE.flipangle == STE.flipangle
        @assert SE.TR == STE.TR
        @assert SE.TM == STE.TM

        new(SE,STE)
    end
end


"""
    WeightedMultiechoContrast([w(TE₁), w(TE₂), ...])

Vector type to store multiecho variable flip angle (VFA) MRI data
"""
struct WeightedMultiechoContrast
    signal
    flipangle
    TR
    TE
    τ
    TM

    # use constructor to check consistency of data
    function WeightedMultiechoContrast(echoList::Vector{WeightedContrast})

        signal = stack((w.signal for w in echoList),dims=1)
        flipangle = first(echoList).flipangle
        TR = first(echoList).TR
        TE = (w.TE for w in echoList)
        τ = first(echoList).τ
        TM = first(echoList).TM

        for w in echoList
            @assert w.flipangle == flipangle "All flip angles must match!"
            @assert w.TR == TR "All TRs must match!"
        end

        new(signal,flipangle,TR,TE,τ,TM)
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