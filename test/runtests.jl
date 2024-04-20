using Antique
using Test
using Suppressor
using Printf
using Markdown
using QuadGK
using Symbolics
using Latexify
using LaTeXStrings
using SpecialFunctions

@testset "Antique.jl" begin
	for model in [:CoulombTwoBody] # Antique.models
		result = @capture_out begin
			include("./$(model).jl")
		end
		open("./result/$(model).log", "w") do io
			println(io, result)
		end
	end
end
