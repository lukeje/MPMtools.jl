using MPMtools.MRIutils
using MPMtools.MRItypes
using MPMtools.MRImaps

# Test MPM pipeline using synthetic data
A      = 1000.0
R1     = 1.0e-3
R2star = 30.0e-3
B1     = 1.0

TR     = 45e-3
TElist = (2:5:40) .* 1e-3
α = MRIutils.optimalDFAangles(TR,R1)
αErnst = MRIutils.ernstangle(TR,R1)

ϵRepeat = MRImaps.dR2star(R1,R2star,[1.0,1.0,1.0],vcat(α...,last(α)),fill(TR,3),fill(TElist,3))
ϵErnst  = MRImaps.dR2star(R1,R2star,[1.0,1.0,1.0],vcat(α...,αErnst),fill(TR,3),fill(TElist,3))

# adding Ernst angle measurement should improve R2star estimate
@test first(ϵErnst) < first(ϵRepeat)

# adding Ernst angle measurement should improve T1w and PDw signal estimate
@test first(last(ϵErnst),2) < first(last(ϵRepeat),2)