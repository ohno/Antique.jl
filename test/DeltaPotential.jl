DP = DeltaPotential(α=1.0, m=1.0, ℏ=1.0)

# <ψᵢ|ψⱼ> = ∫ψ*ψdx = δᵢⱼ

println(raw"""
#### Normalization of $\psi(x)$

```math
\int_{-\infty}^{\infty} \psi^\ast(x) \psi(x) ~\mathrm{d}x = 1
```

```""")

@testset "<ψ|ψ> = ∫ψ*ψdx = 1" begin
  println("  α |   m |   ℏ |        analytical |         numerical ")
  println("--- | --- | --- | ----------------- | ----------------- ")
  for α in [0.1, 1.0, 7.0]
  for m in [0.1, 1.0, 7.0]
  for ℏ in [0.1, 1.0, 7.0]
    DP = DeltaPotential(α=α, m=m, ℏ=ℏ)
    analytical = 1
    numerical  = quadgk(x -> conj(ψ(DP, x)) * ψ(DP, x), -Inf, Inf, maxevals=10^3, order=10)[1]
    acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-5) : isapprox(analytical, numerical, rtol=1e-5)
    @test acceptance
    @printf("%.1f | %.1f | %.1f | %17.12f | %17.12f %s\n", α, m, ℏ, analytical, numerical, acceptance ? "✔" : "✗")
  end
  end
  end
end

println("""```
""")

