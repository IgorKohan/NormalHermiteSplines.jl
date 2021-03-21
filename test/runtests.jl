using NormalHermiteSplines
using DoubleFloats
using Test

@testset "NormalHermiteSplines.jl" begin

include("1D.jl")
include("2D.jl")
include("3D.jl")

end
