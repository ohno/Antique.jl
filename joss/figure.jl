# Rayleigh-Ritz method
using TwoBody
H = Hamiltonian(NonRelativisticKinetic(1, 1), CoulombPotential(-1))
BS = GeometricBasisSet(SimpleGaussianBasis, 0.1, 80.0, 20)
res = solve(H, BS; info=0)
E = res.E
C = res.C
S = res.S

# benchmark
using Antique: Antique
HA = Antique.HydrogenAtom(; Z=1, Eₕ=1.0, a₀=1.0, mₑ=1.0, ℏ=1.0)

# print
for n = 1:length(BS)
  println("\nn = $n\n")
  println("Norm")
  println("  numerical : ", transpose(C[:, n]) * S * C[:, n])
  println("  analytical: ", 1)
  println("Energy")
  println("  numerical : ", E[n])
  println("  analytical: ", Antique.E(HA; n=n))
end

# plot
using CairoMakie
fig = Figure(; size=(840, 600), fontsize=11, backgroundcolor=:transparent)
for i = 1:4
  n = [1, 2, 6, 8][i]
  xₘₐₓ = [6, 18, 120, 200][i]
  yₘₐₓ = 1.06 * maximum(4π * r^2 * abs(Antique.ψ(HA, r, 0, 0; n=n))^2 for r = 0:0.1:xₘₐₓ)
  axis = Axis(fig[div(i - 1, 2)+1, rem(i - 1, 2)+1]; backgroundcolor=:transparent, xlabel=L"$r~/~a_0$", ylabel=L"$4\pi r^2|\psi(r)|^2~ /~{a_0}^{-1}$", xlabelsize=16.5, ylabelsize=16.5, limits=(0, xₘₐₓ, 0, yₘₐₓ))
  lines!(axis, 0 .. 250, r -> 4π * r^2 * abs(res.ψ[n](r))^2; label="Numerical")
  lines!(axis, 0 .. 250, r -> 4π * r^2 * abs(Antique.ψ(HA, r, 0, 0; n=n))^2; label="Analytical", color=:black, linestyle=:dash)
  axislegend(axis, "n = $n"; position=:rt, framevisible=false)
end
file = String(@__FILE__)
save(replace(file, r"(?<path>.*)\.jl" => s"\g<path>.pdf"), fig)
fig
