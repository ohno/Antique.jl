io = open("./result/DeltaPotential.md", "w")


# <ψᵢ|ψⱼ> = ∫ψ(x)*ψ(x)dx = δᵢⱼ


println(io, raw"""
#### Normalization of $\psi(x)$

```math
\int_{-\infty}^{\infty} \psi^\ast(x) \psi(x) ~\mathrm{d}x = 1
```

```""")

@testset "DP: <ψ|ψ> = ∫ψ(x)*ψ(x)dx = 1" begin
  println(io, "  α |   m |   ℏ |     analytical |      numerical ")
  println(io, "--- | --- | --- | -------------- | -------------- ")
  for α in [0.1, 1.0, 7.0]
  for m in [0.1, 1.0, 7.0]
  for ℏ in [0.1, 1.0, 7.0]
    DP = DeltaPotential(α=α, m=m, ℏ=ℏ)
    analytical = 1
    numerical  = quadgk(x -> conj(ψ(DP, x)) * ψ(DP, x), -Inf, Inf, maxevals=10^3, order=10)[1]
    acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-5) : isapprox(analytical, numerical, rtol=1e-5)
    @test acceptance
    @printf(io, "%.1f | %.1f | %.1f | %14.9f | %14.9f %s\n", α, m, ℏ, analytical, numerical, acceptance ? "✔" : "✗")
  end
  end
  end
end

println(io, """```\n""")


# <ψₙ|H|ψₙ>  = ∫ψₙ*Tψₙdx = Eₙ


println(io, raw"""
#### Eigenvalues

```math
\begin{aligned}
  E_n
  &= \int_{-\infty}^{\infty} \psi^\ast(x) \hat{H} \psi(x) ~\mathrm{d}x \\
  &= \int_{-\infty}^{\infty} \psi^\ast(x) \left[ \hat{V} + \hat{T} \right] \psi(x) ~\mathrm{d}x \\
  &= \int_{-\infty}^{\infty} \psi^\ast(x) \left[ -\alpha\delta(x) - \frac{\hbar^2}{2m} \frac{\mathrm{d}^{2}}{\mathrm{d} x^{2}} \right] \psi(x) ~\mathrm{d}x \\
  &= \int_{-\infty}^{\infty} \psi^\ast(x) \left[ -\alpha\delta(x) - \frac{\hbar^2}{2m} (\kappa^2 -2\kappa \delta(x))\right]\psi(x) ~\mathrm{d}x \\
\end{aligned}
```
where the $\kappa=m\alpha/\hbar^2$ and the integration with the delta function yeild a term proportional to the wave function at the origin $|\psi(0)|^2$.
```""")

@testset "DP: <ψₙ|H|ψₙ> = ∫ψₙ*Hψₙdx = Eₙ" begin
  function ∫ψHψdx(DP)
    κ = DP.m * DP.α/ DP.ℏ^2
    T = -DP.ℏ^2 / (2 * DP.m) *( κ^2 - 2κ * ψ(DP,0)^2 )
    V = -DP.α * ψ(DP, 0)^2
    return T + V
  end
  println(io, "  α |   m |   ℏ |     analytical |      numerical ")
  println(io, "--- | --- | --- | -------------- | -------------- ")
  for α in [0.1, 1.0, 7.0]
  for m in [0.1, 1.0, 7.0]
  for ℏ in [0.1, 1.0, 7.0]
    DP = DeltaPotential(α=α, m=m, ℏ=ℏ)
    analytical = E(DP)
    numerical = ∫ψHψdx(DP)
    acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-4) : isapprox(analytical, numerical, rtol=1e-4)
    @test acceptance
    @printf(io, "%.1f | %.1f | %.1f | %14.9f | %14.9f %s\n", α, m, ℏ,analytical, numerical, acceptance ? "✔" : "✗") 
  end
  end
  end
end

println(io, """```\n""")


close(io)