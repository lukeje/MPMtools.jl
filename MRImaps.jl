module MRImaps

using LinearAlgebra

"""
    WeightedContrast(signal, flipangle, TR[, B1])

Composite type to store variable flip angle (VFA) MRI data
"""
struct WeightedContrast
    signal
    flipangle::Number
    TR::Number
    B1
    τ

    # define explicit inner constructor as τ is computed from the input
    # and we need to check some parameter bounds
    function WeightedContrast(signal, flipangle, TR, B1)

        @assert (TR > 0) "TR must be greater than zero!"

        # populate τ field (half angle tangent transform) using provided B1 and flipangle
        τ = 2tan(B1 * flipangle / 2)

        new(signal, flipangle, TR, B1, τ)
    end
end

# outer constructor: use B1 = 1.0 if B1 not provided
WeightedContrast(signal, flipangle, TR) = WeightedContrast(signal, flipangle, TR, 1.0)

"""
    WeightedMultiechoContrast([w(TE₁), w(TE₂), ...], [TE₁, TE₂, ...])

Composite type to store multiecho variable flip angle (VFA) MRI data
"""
struct WeightedMultiechoContrast
    contrastList::Vector{WeightedContrast}
    TElist::Vector{Number}

    function WeightedMultiechoContrast(contrastList, TElist)
        @assert (length(contrastList) == length(TElist)) "Number of contrasts must match provided number of echoes!"
        new(contrastList,TElist)
    end
end

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

"""
    calculateR2star(weighted_data)

R2* estimation using an implementation of the ESTATICS model (Weiskopf2014)

# Arguments
array of WeightedMultiechoContrast (one per contrast)
- Voxels must correspond between the weightings (i.e. the images should
    have been resliced to the same space), but the sampled TEs may be
    different.
- nTEs must be at least 2 for each weighting.
- Because log(0) is ill-defined, zero values in any voxel will result
    in NaN output for that voxel. To avoid potentially biasing the data,
    we do not modify the input in any way to avoid this, and leave it to
    the user to decide how to handle this case, e.g. by removing the
    corresponding voxels from the input data or replacing zeroes with a
    small positive number.

# Outputs
R2star (nVoxelsX x nVoxelsY x ...): the voxelwise-estimated common R2* of the weightings.

# Examples
ESTATICS estimate of common R2* of PD, T1, and MT-weighted data:
```juliarepl
TBD
```

# References:
Weiskopf et al. Front. Neurosci. (2014), "Estimating the apparent transverse relaxation time (R2*) from images with different contrasts
(ESTATICS) reduces motion artifacts", [doi:10.3389/fnins.2014.00278](https://doi.org/10.3389/fnins.2014.00278)
"""
function calculateR2star(weighted_dataList::Vector{WeightedMultiechoContrast})
        
    dims = size(weighted_dataList[1].contrastList[1].signal)
    nVoxels = prod(dims)
    nWeighted = length(weighted_dataList)
    
    # Build regression arrays
    # Build design matrix
    D = Array{Float64}(undef, 0, nWeighted+1)
    for w = 1:nWeighted
        d = zeros(length(weighted_dataList[w].TElist), nWeighted+1)
        d[:,1] = -weighted_dataList[w].TElist
        d[:,w+1] .= 1
        D = vcat(D, d)
    end
    
    # Build response variable vector
    y = Array{Float64}(undef, 0, 1)
    for w = 1:nWeighted
        
        nTEs = length(weighted_dataList[w].TElist)
        @assert (nTEs > 1) "each weighting must have more than one TE"
        
        localDims=size(weighted_dataList[w].contrastList[1].signal)
        @assert (prod(localDims) == nVoxels) "all input data must have the same number of voxels"
        
        rData = zeros(nVoxels, nTEs)
        for t = 1:nTEs
            rData[:, t] .= weighted_dataList[w].contrastList[t].signal
        end
        
        # log(0) is not defined, so warn the user about zeroes in their data 
        # for methods involving a log transform.
        if any(rData .== 0)
            @warn """Zero values detected in some voxels in the input data. This will cause estimation to fail in these voxels due to the log 
                transform. If these voxels are background voxels, consider removing them from the input data matrices. 
                Zero values which occur only at high TE in voxels of interest could be replaced with a small positive number, e.g. eps()
                (if the data magnitudes are ≈ 1) or 1 if the data are integer-valued. Note: Care must be taken when replacing 
                values, as this could bias the R2* estimation."""
        end
        
        y = vcat(y, transpose(rData))
    end
    
    # Estimate R2*
    β = (transpose(D) * D) \ (transpose(D) * log.(y))
    β[2:end,:] = exp.(β[2:end,:])
    
    # Output
    # extra unity in reshape argument avoids problems if size(dims)==2
    return R2s = reshape(β[1,:],dims)
end

end