using Antique
using ForwardDiff
using HCubature
using Latexify
using LaTeXStrings
using Markdown
using SpecialFunctions
using Symbolics
using Test
using Printf
using QuadGK
using Zygote

@testset verbose = true "Antique.jl" begin
	for model in Antique.models # [:InfinitePotentialWell3D, :SphericalOscillator, :HydrogenAtom, :CoulombTwoBody]
		include("./$(model).jl")
	end
end
