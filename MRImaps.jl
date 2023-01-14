module MRImaps

"""
    WeightedContrast(signal, flipangle, TR, TE[, B1])

Composite type to store variable flip angle (VFA) MRI data
"""
struct WeightedContrast
    signal
    flipangle::Number
    TR::Number
    TE::Union{Number, Missing}
    B1
    τ

    # define explicit inner constructor as τ is computed from the input
    # and we need to check some parameter bounds
    function WeightedContrast(signal, flipangle, TR, TE, B1)

        @assert (TR > 0) "TR must be greater than zero!"

        # populate τ field (half angle tangent transform) using provided B1 and flipangle
        τ = 2tan(B1 * flipangle / 2)

        new(signal, flipangle, TR, TE, B1, τ)
    end
end

# outer constructor: use B1 = 1.0 if B1 not provided
WeightedContrast(signal, flipangle, TR, TE) = WeightedContrast(signal, flipangle, TR, TE, 1.0)

"""
    calculateA(PDw::WeightedContrast, T1w::WeightedContrast)

Calculate A (an unnormalised PD map) from PDw and T1w signals
 
# Arguments
PDw and T1w are both of type WeightedContrast

# Output
A (in arbitrary units)

# Examples
```juliarepl
PDw = WeightedContrast(signal_pdw, fa_pdw, tr_pdw, b1map)
T1w = WeightedContrast(signal_t1w, fa_t1w, tr_t1w, b1map)
A = calculateA(PDw, T1w)
```

# References
Helms et al. Magn. Reson. Med. (2008), "Quantitative FLASH MRI at 3T using a rational approximation of the Ernst equation", [doi:10.1002/mrm.21542](https://doi.org/10.1002/mrm.21542)

Edwards et al.  Magn. Reson. Mater. Phy. (2021), "Rational approximation of the Ernst equation for dual angle R1 mapping revisited: beyond the small flip-angle assumption" 
in Book of Abstracts ESMRMB 2021, [doi:10.1007/s10334-021-00947-8](https://doi.org/10.1007/s10334-021-00947-8)
"""
function calculateA(PDw::WeightedContrast, T1w::WeightedContrast)
    
    @assert all(size(PDw.signal)==size(T1w.signal)) "PDw.signal and T1w.signal must be the same size!"
    
    if (PDw.signal ≠ 0) & (T1w.signal ≠ 0) & (PDw.τ ≠ 0) & (T1w.τ ≠ 0)
        return PDw.signal * T1w.signal * ( T1w.TR * PDw.τ / T1w.τ - PDw.TR * T1w.τ / PDw.τ )  / ( PDw.signal * T1w.TR * PDw.τ - T1w.signal * PDw.TR * T1w.τ )
    else # cannot get sensible answer
        return missing
    end
end

"""
    calculateR1(PDw::WeightedContrast, T1w::WeightedContrast)

Calculate R1 map from PDw and T1w signals

# Arguments
PDw and T1w are both of type WeightedContrast

# Output
R1 (in reciprocal of TR units)

# Examples
```juliarepl
PDw = WeightedContrast(signal_pdw, fa_pdw, tr_pdw, b1map)
T1w = WeightedContrast(signal_t1w, fa_t1w, tr_t1w, b1map)
R1 = calculateR1(PDw, T1w)
```

# References
Helms et al. Magn. Reson. Med. (2008), "Quantitative FLASH MRI at 3T using a rational approximation of the Ernst equation", [doi:10.1002/mrm.21542](https://doi.org/10.1002/mrm.21542)

Edwards et al.  Magn. Reson. Mater. Phy. (2021), "Rational approximation of the Ernst equation for dual angle R1 mapping revisited: beyond the small flip-angle assumption" 
in Book of Abstracts ESMRMB 2021, [doi:10.1007/s10334-021-00947-8](https://doi.org/10.1007/s10334-021-00947-8)
"""
function calculateR1(PDw::WeightedContrast, T1w::WeightedContrast)

    @assert all(size(PDw.signal)==size(T1w.signal)) "PDw.signal and T1w.signal must be the same size!"
        
    if (PDw.signal ≠ 0) & (T1w.signal ≠ 0) & (PDw.τ ≠ 0) & (T1w.τ ≠ 0)
        return  0.5 * ( PDw.signal * PDw.τ / PDw.TR - T1w.signal * T1w.τ / T1w.TR ) / ( T1w.signal / T1w.τ - PDw.signal / PDw.τ )
    else # cannot get sensible answer
        return missing
    end
end

end