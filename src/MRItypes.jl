module MRItypes

export WeightedContrast, WeightedMultiechoContrast

"""
    WeightedContrast(signal, flipangle, TR, B1[, TE])

Composite type to store variable flip angle (VFA) MRI data
"""
struct WeightedContrast
    signal
    flipangle::Number
    TR::Number
    B1
    TE::Union{Number,Missing}
    τ

    # define explicit inner constructor as τ is computed from the input
    # and we need to check some parameter bounds
    function WeightedContrast(signal, flipangle, TR, B1, TE)

        @assert (TR > 0) "TR must be greater than zero!"
        @assert all((B1 .> 0) .| ismissing(B1)) "B1 must be greater than zero or missing!"
        @assert (TE ≥ 0) | ismissing(TE) "TE must be greater than or equal to zero or missing!"

        # populate τ field (half angle tangent transform) using provided B1 and flipangle
        τ = 2tan(B1 * flipangle / 2)

        new(signal, flipangle, TR, B1, TE, τ)
    end
end

# outer constructor: TE is missing if not provided
WeightedContrast(signal, flipangle, TR, B1) = WeightedContrast(signal, flipangle, TR, B1, missing)

"""
    WeightedMultiechoContrast([w(TE₁), w(TE₂), ...])

Vector type to store multiecho variable flip angle (VFA) MRI data
"""
struct WeightedMultiechoContrast
    echoList::Vector{WeightedContrast}
    # Really just an alias of the vector type; use constructor to check consistency of data
    function WeightedMultiechoContrast(echoList::Vector{WeightedContrast})
        for w in echoList
            @assert w.flipangle == echoList[1].flipangle "All flip angles must match!"
            @assert w.TR == echoList[1].TR "All TRs must match!"
            @assert w.TE > 0 "All TEs must be greater than zero!"
            @assert size(w.signal) == size(echoList[1].signal) "The dimensions of the data must match!"
        end
        @assert (length(unique([c.TE for c in echoList])) > 1) "There must be more than one unique TE!"
        
        new(echoList)
    end
end

end