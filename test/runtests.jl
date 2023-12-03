using AnalyticalSolutions
using Test
using Suppressor

@testset "AnalyticalSolutions.jl" begin
	@suppress_out begin
		include("./InfinitePotentialWell.jl")
		include("./HarmonicOscillator.jl")
		include("./MorsePotential.jl")
		include("./HydrogenAtom.jl")
	end
end
