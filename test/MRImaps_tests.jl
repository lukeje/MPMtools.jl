using MPMtools.MRIutils
using MPMtools.MRItypes
using MPMtools.MRImaps

# Test MPM pipeline using synthetic data
A      = 1000.0
R1     = 1.0
R2star = 30.0
B1     = 1.0

TR     = 50e-3
TElist = (2:5:40) .* 1e-3
faPDw  = B1*6
faT1w  = B1*18

PDw_A = MRIutils.ernstd(faPDw, TR, R1, PD=A)
T1w_A = MRIutils.ernstd(faT1w, TR, R1, PD=A)

PDw = MRIutils.exponentialDecay(PDw_A, R2star, TElist)
T1w = MRIutils.exponentialDecay(T1w_A, R2star, TElist)

# log linear ordinary least squares
R2star_OLS_est, (PDw0, T1w0) = MRImaps.calculateR2star([PDw, T1w])
@test R2star_OLS_est[] ≈ R2star

A_est  = MRImaps.calculateA(PDw0, T1w0)
@test A_est[] ≈ A

R1_est = MRImaps.calculateR1(PDw0, T1w0)
@test R1_est[] ≈ R1 rtol=1e-3

(A_est,T1_est) = MRImaps.calculateT1([PDw0, T1w0])

@test A_est[] ≈ A
@test T1_est[] ≈ (1/R1) rtol=1e-3

# iterative log linear weighted least squares
R2star_WLS_est, _ = MRImaps.calculateR2star([PDw, T1w],niter=3)

@test R2star_WLS_est[] ≈ R2star
