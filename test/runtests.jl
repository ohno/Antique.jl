using Antique
using Test
using Suppressor

@testset "Antique.jl" begin
	@suppress_out begin
		for model in Antique.models
			include("./$(model).jl")
		end
	end
end
