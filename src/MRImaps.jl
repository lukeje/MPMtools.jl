module MRImaps

using LinearAlgebra
using ..MRItypes

"""
    calculateA(PDw::WeightedContrast, T1w::WeightedContrast)

Calculate A (an unnormalised PD map) from PDw and T1w signals
 
# Arguments
PDw and T1w are both of type WeightedContrast

# Output
A (in arbitrary units)

# Examples
```juliarepl
PDw = WeightedContrast(signal_pdw, fa_pdw, tr_pdw)
T1w = WeightedContrast(signal_t1w, fa_t1w, tr_t1w)
A = calculateA(PDw, T1w)
```

# References
- Helms et al. Magn. Reson. Med. (2008), "Quantitative FLASH MRI at 3T using a rational approximation of the Ernst equation", [doi:10.1002/mrm.21542](https://doi.org/10.1002/mrm.21542)
- Edwards et al.  Magn. Reson. Mater. Phy. (2021), "Rational approximation of the Ernst equation for dual angle R1 mapping revisited: beyond the small flip-angle assumption" in Book of Abstracts ESMRMB 2021, [doi:10.1007/s10334-021-00947-8](https://doi.org/10.1007/s10334-021-00947-8)
"""
function calculateA(PDw::WeightedContrast, T1w::WeightedContrast)
    
    @assert all(size(PDw.signal)==size(T1w.signal)) "PDw.signal and T1w.signal must be the same size!"
    @assert (PDw.TR ≠ T1w.TR) | (PDw.τ ≠ T1w.τ) "PDw and T1w data must differ in either TR or flip angle!"
    
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
PDw = WeightedContrast(signal_pdw, fa_pdw, tr_pdw)
T1w = WeightedContrast(signal_t1w, fa_t1w, tr_t1w)
R1 = calculateR1(PDw, T1w)
```

# References
- Helms et al. Magn. Reson. Med. (2008), "Quantitative FLASH MRI at 3T using a rational approximation of the Ernst equation", [doi:10.1002/mrm.21542](https://doi.org/10.1002/mrm.21542)
- Edwards et al.  Magn. Reson. Mater. Phy. (2021), "Rational approximation of the Ernst equation for dual angle R1 mapping revisited: beyond the small flip-angle assumption" in Book of Abstracts ESMRMB 2021, [doi:10.1007/s10334-021-00947-8](https://doi.org/10.1007/s10334-021-00947-8)
"""
function calculateR1(PDw::WeightedContrast, T1w::WeightedContrast)

    @assert all(size(PDw.signal)==size(T1w.signal)) "PDw.signal and T1w.signal must be the same size!"
    @assert (PDw.TR ≠ T1w.TR) | (PDw.τ ≠ T1w.τ) "PDw and T1w data must differ in either TR or flip angle!"
        
    if (PDw.signal ≠ 0) & (T1w.signal ≠ 0) & (PDw.τ ≠ 0) & (T1w.τ ≠ 0)
        return  0.5 * ( PDw.signal * PDw.τ / PDw.TR - T1w.signal * T1w.τ / T1w.TR ) / ( T1w.signal / T1w.τ - PDw.signal / PDw.τ )
    else # cannot get sensible answer
        return missing
    end
end


"""
    calculateT1([PDw::WeightedContrast, T1w::WeightedContrast, ...])

Calculate T1 map from any number of weighted signals

# Arguments
PDw and T1w are both of type WeightedContrast

# Output
T1 (in TR units)

# Examples
```juliarepl
PDw = WeightedContrast(signal_pdw, fa_pdw, tr_pdw)
T1w = WeightedContrast(signal_t1w, fa_t1w, tr_t1w)
(A,T1) = calculateR1([PDw, T1w])
```

# References
- Helms et al. Magn. Reson. Med. (2011), "Identification of signal bias in the variable flip angle method by linear display of the algebraic ernst equation", [doi:10.1002/mrm.22849](https://doi.org/10.1002/mrm.22849)
"""
function calculateT1(weighted::Vector{WeightedContrast})

    @assert all(size(w.signal)==size(first(weighted).signal) for w in weighted) "All weighted data must be the same size!"
    
    S  = [w.signal for w in weighted]
    τ  = [w.τ      for w in weighted]
    TR = [w.TR     for w in weighted]

    y = S./τ
    D = hcat(ones(length(weighted)), -S.*τ./(2TR))

    (A,T1) = D\y

end


"""
    calculateR2star(weighted_data)

R2* estimation using an implementation of the ESTATICS model (Weiskopf2014)

# Arguments
array of WeightedMultiechoContrast (one per contrast)
- Voxels must correspond between the weightings (i.e. the images should have been resliced to the same space), but the sampled TEs may be
    different.
- Because log(0) is ill-defined, zero values in any voxel will result in NaN output for that voxel. To avoid potentially biasing the data,
    we do not modify the input in any way to avoid this, and leave it to the user to decide how to handle this case, e.g. by removing the
    corresponding voxels from the input data or replacing zeroes with a small positive number.

# Outputs
R2star: the estimated common R2* of the weightings.

# Examples
ESTATICS estimate of common R2* of PD, T1, and MT-weighted data:
```juliarepl
TBD
```

# References:
- Weiskopf et al. Front. Neurosci. (2014), "Estimating the apparent transverse relaxation time (R2*) from images with different contrasts (ESTATICS) reduces motion artifacts", [doi:10.3389/fnins.2014.00278](https://doi.org/10.3389/fnins.2014.00278)
"""
function calculateR2star(weighted_dataList::Vector{WeightedMultiechoContrast})

    T = Float64
        
    nWeighted = length(weighted_dataList)
    
    # Build design matrix D and response variable y
    D = Matrix{T}(undef, 0, nWeighted+1)
    y = Vector{T}(undef, 0)
    for wIdx = eachindex(weighted_dataList)
        w = weighted_dataList[wIdx]
        nTEs = length(w.TE)

        d = hcat(-T.(w.TE), zeros(T, nTEs, wIdx-1), ones(T, nTEs, 1), zeros(T, nTEs, nWeighted-wIdx))
        D = vcat(D, d)
        
        append!(y, T.(w.signal))
    end

    # log(0) is not defined, so warn the user about zeroes in their data 
    # for methods involving a log transform.
    if any(y .== 0)
        @warn """Zero values detected in the input data. This will cause estimation to fail in these voxels due to the log 
            transform. If these voxels are background voxels, consider removing them from the input data matrices. 
            Zero values which occur only at high TE in voxels of interest could be replaced with a small positive number, e.g. eps()
            (if the data magnitudes are ≈ 1) or 1 if the data are integer-valued. Note: Care must be taken when replacing 
            values, as this could bias the R2* estimation."""
    end
    
    # Estimate R2*
    β = (transpose(D) * D) \ (transpose(D) * log.(y))

    # Output
    R2star = β[1]

    extrapolated = Vector{WeightedContrast}(undef,0)
    for wIdx = eachindex(weighted_dataList)
        w = weighted_dataList[wIdx]
        push!(extrapolated, WeightedContrast(exp(β[wIdx+1]), w.flipangle, w.TR, zero(w.TE[1])))
    end

    return R2star, extrapolated
end


end
