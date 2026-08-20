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

    @testset "Base.string" begin
        @test string(HarmonicOscillator(k = 2.0, m = 3.0, hbar = 4.0)) == "Antique.HarmonicOscillator(k=2.0, m=3.0, hbar=4.0)"
        for model_name in Antique.models
            model = getproperty(Antique, model_name)()
            representation = string(model)
            @test representation == sprint(show, model)
            @test startswith(representation, "Antique.$(model_name)(")
            @test which(Base.string, (typeof(model),)).module === Base
        end
    end

    for model in Antique.models # [:InfinitePotentialWell3D, :SphericalOscillator, :HydrogenAtom, :CoulombTwoBody]
        include("./$(model).jl")
    end
end
