using MPMtools
using Test

@testset "MRImaps tests" begin
    include("MRImaps_tests.jl")
    include("MRImaps_errormaps_tests.jl")
end