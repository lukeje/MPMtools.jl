module MRImaps

using LinearAlgebra
using ..MRItypes
using ..MRIutils
using Statistics: std, mean
using Combinatorics: with_replacement_combinations

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
    @assert (PDw.TR ≠ T1w.TR) || (PDw.τ ≠ T1w.τ) "PDw and T1w data must differ in either TR or flip angle!"

    return @. PDw.signal * T1w.signal * ( T1w.TR * PDw.τ / T1w.τ - PDw.TR * T1w.τ / PDw.τ )  / ( PDw.signal * T1w.TR * PDw.τ - T1w.signal * PDw.TR * T1w.τ )
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
    @assert (PDw.TR ≠ T1w.TR) || (PDw.τ ≠ T1w.τ) "PDw and T1w data must differ in either TR or flip angle!"

    return @. 0.5 * ( PDw.signal * PDw.τ / PDw.TR - T1w.signal * T1w.τ / T1w.TR ) / ( T1w.signal / T1w.τ - PDw.signal / PDw.τ )
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
    
    #TODO: include case where τ varies per voxel
    S  = stack((w.signal for w in weighted), dims=1)
    τ  = stack((w.τ      for w in weighted), dims=1)
    TR = stack((w.TR     for w in weighted), dims=1)

    A  = similar(first(weighted).signal)
    T1 = similar(A)
    for i in CartesianIndices(A)
        (A[i],T1[i]) = _calculateT1(S[:,i],τ,TR)
    end

    return A,T1

end
function _calculateT1(S,τ,TR)
    y = S./τ
    D = hcat(ones(eltype(S),length(S)), -S.*τ./(2TR))

    (A,T1) = D\y
end

"""
    calculateR2star(weighted_data)

R2* estimation using an implementation of the ESTATICS model (Weiskopf2014)

# Arguments
- array of WeightedMultiechoContrast (one per contrast)
    - Voxels must correspond between the weightings (i.e. the images should have been resliced to the same space), but the sampled TEs may be
      different.
    - Because log(0) is ill-defined, zero values in any voxel will result in NaN output for that voxel. To avoid potentially biasing the data,
      we do not modify the input in any way to avoid this, and leave it to the user to decide how to handle this case, e.g. by removing the
      corresponding voxels from the input data or replacing zeroes with a small positive number.
- niter: number of iterations if weighted least squares is to be used. Default is 0 (ordinary least squares)

# Outputs
R2star: the estimated common R2* of the weightings.

# Examples
ESTATICS estimate of common R2* of PD, T1, and MT-weighted data:
```juliarepl
TBD
```

# References:
- Weiskopf et al. Front. Neurosci. (2014), "Estimating the apparent transverse relaxation time (R2*) from images with different contrasts (ESTATICS) reduces motion artifacts", [doi:10.3389/fnins.2014.00278](https://doi.org/10.3389/fnins.2014.00278)
- Edwards et al. Proc. Int. Soc. Magn. Reson. Med. (2022), "Robust and efficient R2* estimation in human brain using log-linear weighted least squares"
"""
function calculateR2star(weighted_dataList::Vector{WeightedMultiechoContrast}; niter=0)

    T = Float64

    nWeighted = length(weighted_dataList)
    nVoxels = size(first(weighted_dataList).signal)[2:end]

    # Build design matrix D and response variable y
    D = Matrix{T}(undef, 0, nWeighted+1)
    for (i,w) in enumerate(weighted_dataList)
        nTEs = length(w.TE)

        d = hcat(-T.(w.TE), zeros(T, nTEs, i-1), ones(T, nTEs, 1), zeros(T, nTEs, nWeighted-i))
        D = vcat(D, d)

        @assert size(w.signal)[2:end] == nVoxels "all weighted data must have the same number of voxels"
    end

    y = reduce(vcat, (T.(w.signal) for w in weighted_dataList))

    # log(0) is not defined, so warn the user about zeroes in their data
    if any(y .== 0)
        @warn """Zero values detected in the input data. This will cause estimation to fail in these voxels due to the log 
            transform. If these voxels are background voxels, consider removing them from the input data matrices. 
            Zero values which occur only at high TE in voxels of interest could be replaced with a small positive number, e.g. eps()
            (if the data magnitudes are ≈ 1) or 1 if the data are integer-valued. Note: Care must be taken when replacing 
            values, as this could bias the R2* estimation."""
    end

    # Estimate R2* using iterative weighted least squares
    R2star = Array{T}(undef,(isempty(nVoxels) ? 1 : nVoxels)...)
    extrapolated = Array{T}(undef,nWeighted,(isempty(nVoxels) ? 1 : nVoxels)...)
    for v in CartesianIndices(size(y)[begin+1:end])
        yloc = view(y,:,v)

        if any(isnan.(yloc))
            R2star[v] = T(NaN)
            extrapolated[:,v] .= T(NaN)
        else
            S = ones(T,length(yloc)) # run at least one iteration with OLS
            logy = log.(yloc)
            for iter in 0:niter # zeroth iteration is OLS
                S2 = S.^2
                
                # return earlier iteration (usually OLS result) if all signals are zero
                # in current iteration
                sum(S2) == zero(T) && break

                W = diagm(S2)

                global β = (transpose(D) * W * D) \ (transpose(D) * W * logy)
                S = exp.(D*β)
            end

            # Output
            R2star[v] = β[begin]
            extrapolated[:,v] = exp.(β[begin+1:end])
        end
    end

    return R2star, [WeightedContrast(b, w.flipangle, w.TR, zero(first(w.TE))) for (w,b) in zip(weighted_dataList,eachslice(extrapolated,dims=1))]
end


"""
    dR2star([SPD,ST1,...], [dSPD,dST1,...])
Calculate propagation of uncertainty for T1 and A map.

In: ...

Out: ...

# References
- https://en.wikipedia.org/wiki/Propagation_of_uncertainty
- ...
    
"""
function dR2star(weighted_dataList, dS, S0)
    
    T = Float64
        
    nweighted = length(weighted_dataList)
    
    # Build design matrix D and response variable y
    D   = Matrix{T}(undef, 0, nweighted+1)
    y   = Vector{T}(undef, 0)
    dSy = Vector{T}(undef, 0)
    for wIdx = eachindex(weighted_dataList)
        w = weighted_dataList[wIdx]
        nTEs = length(w.TE)

        d = hcat(-T.(w.TE), zeros(T, nTEs, wIdx-1), ones(T, nTEs, 1), zeros(T, nTEs, nweighted-wIdx))
        D = vcat(D, d)
        
        append!(y,   T.(w.signal))
        append!(dSy, T.(fill(dS[wIdx],nTEs)))
    end

    d = zeros(T,nweighted + 1)
    for (n,yloc) = pairs(y)
        logy′ = zero(y)
        logy′[n] = one(yloc)/yloc

        # use OLS solution if trace of weight matrix zero
        # to match R2star calculation
        y2 = y.^2
        sum(y2) == 0 ? W = I : W = diagm(y2) ./ sum(y2)

        # ignore dependence of weights on signal for now
        σ = (transpose(D) * W * D) \ (transpose(D) * W * logy′)

        d .+= (σ .* dSy[n]).^2
    end

    return sqrt(d[1]), sqrt.(S0.signal^2 * d for (S0,d) in zip(S0,d[begin+1:end]))
end

function dR2star(weighted_dataList, dS)
    _,S0 = calculateR2star(weighted_dataList; niter=0)
    dR2star(weighted_dataList,dS,S0)
end

function dR2star(R₁,R2star,dS,α,TR,TE)
    S0 = (MRIutils.ernst(α,TR,R₁) for (α,TR) in zip(α,TR))
    dR2star([MRIutils.exponentialDecay(S0,R2star,TE) for (S0,TE) in zip(S0,TE)], dS, S0)
end


"""
    calculateSESTEB1(W, nominalR1, nse=5, nambiguousangles=2)
Calculate B1 estimates for data acquired using a spin echo (SE)/stimulated echo (STE) protocol.
Expects as input a vector of SE/STE contrast pairs each acquired with a different nominal flip 
angle. A typical R1 value `nominalR1` is used in the correction of R1 decay during the mixing 
time. The first `nse` of these pairs with the highest SE intensities will be used for the 
calculation. As the data is typically magnitude only, the actual angles are ambiguous. The code
will check for up to `nwraps` wraps in the data. Setting this value too high will massively
slow down the computation and lead to very high background values, so `nwraps` = 1 or 2 is a
good choice, depending on the expected range of inhomogeneity.

"""
function calculateSESTEB1(W::Vector{SESTEcontrast}, nominalR1, nse=5, nwraps=2)
    nominalfa = [w.SE.flipangle for w in W]
    wrappedfa = stack((@. real(acos(complex(exp(w.STE.TM*nominalR1) * w.STE.signal / w.SE.signal))) for w in W), dims=1) ./ nominalfa

    B1 = similar(first(W).SE.signal)
    p = floor(-nwraps/2):floor(nwraps/2)
    for c in CartesianIndices(first(W).SE.signal)
        σ = Inf64
        ihise = sortperm([w.SE.signal[c] for w in W],rev=true)[begin:nse]
        y = view(wrappedfa,ihise,c)
        for p in with_replacement_combinations(p,nse)
            y′ = @. abs(y + p*π/nominalfa[ihise])
            if std(y′) < σ
                B1[c] = mean(y′)
                σ = std(y′)
            end
        end
    end
    return B1
end


end
