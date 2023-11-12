using Antiq
using Test

@testset "Antiq.jl" begin
	include("./InfinitePotentialWell.jl")
	include("./HarmonicOscillator.jl")
	include("./MorsePotential.jl")
	include("./HydrogenAtom.jl")
end
