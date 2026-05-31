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
	@testset "Field alias compatibility" begin
		ha_ascii = HydrogenAtom(m_e=2.0, a_0=3.0, E_h=4.0, hbar=5.0)
		@test ha_ascii.mₑ == 2.0
		@test ha_ascii.a₀ == 3.0
		@test ha_ascii.Eₕ == 4.0
		@test ha_ascii.ℏ == 5.0
		ha_unicode = HydrogenAtom(mₑ=6.0, a₀=7.0, Eₕ=8.0, ℏ=9.0)
		@test ha_unicode.m_e == 6.0
		@test ha_unicode.a_0 == 7.0
		@test ha_unicode.E_h == 8.0
		@test ha_unicode.hbar == 9.0

		harm = HarmonicOscillator(hbar=2.0)
		@test harm.ℏ == 2.0

		ctb = CoulombTwoBody(z_1=-2, z_2=3, m_1=1.1, m_2=2.2, m_e=3.3, a_0=4.4, E_h=5.5, hbar=6.6)
		@test ctb.z₁ == -2
		@test ctb.z₂ == 3
		@test ctb.m₁ == 1.1
		@test ctb.m₂ == 2.2
		@test ctb.mₑ == 3.3
		@test ctb.a₀ == 4.4
		@test ctb.Eₕ == 5.5
		@test ctb.ℏ == 6.6

		dp = DeltaPotential(alpha=2.5, hbar=3.5)
		@test dp.α == 2.5
		@test dp.ℏ == 3.5

		ipw = InfinitePotentialWell(hbar=2.5)
		@test ipw.ℏ == 2.5

		ipw3d = InfinitePotentialWell3D(hbar=2.5)
		@test ipw3d.ℏ == 2.5

		mp = MorsePotential(r_e=1.1, D_e=2.2, mu=3.3, hbar=4.4)
		@test mp.rₑ == 1.1
		@test mp.Dₑ == 2.2
		@test mp.μ == 3.3
		@test mp.ℏ == 4.4

		pt = PoschlTeller(lambda=4, x_0=2.5, hbar=3.5)
		@test pt.λ == 4
		@test pt.x₀ == 2.5
		@test pt.ℏ == 3.5

		rr = RigidRotor(m_1=1.2, m_2=2.3, hbar=3.4)
		@test rr.m₁ == 1.2
		@test rr.m₂ == 2.3
		@test rr.ℏ == 3.4

		so = SphericalOscillator(mu=2.2, hbar=3.3)
		@test so.μ == 2.2
		@test so.ℏ == 3.3
	end

	for model in Antique.models # [:InfinitePotentialWell3D, :SphericalOscillator, :HydrogenAtom, :CoulombTwoBody]
		include("./$(model).jl")
	end
end
