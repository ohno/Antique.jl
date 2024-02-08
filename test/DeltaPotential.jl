DP = DeltaPotential(α=1.0,m=1.0, ℏ=1.0)

# <ψᵢ|ψⱼ> = ∫ψ*ψdx = δᵢⱼ

println(raw"""
#### Normalization & Orthogonality of $\psi(x)$

```math
\int_{-Inf}^{Inf} \psi^\ast(x) \psi(x) ~\mathrm{d}x = 1
```

```""")

@testset "<ψ|ψ> = ∫ψ*ψdx = 1" begin
  println("        analytical |         numerical ")
  println(" ----------------- | ----------------- ")
  # for L in [0.1, 0.5, 1.0, 7.0]
  # for m in [0.1, 0.5, 1.0, 7.0]
  # for ℏ in [0.1, 0.5, 1.0, 7.0]
  #  for i in 1:10
  #  for j in 1:10
    analytical = 1
    numerical  = quadgk(x -> conj(ψ(DP, x)) * ψ(DP, x), -Inf, Inf, maxevals=10^3)[1]
    acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-5) : isapprox(analytical, numerical, rtol=1e-5)
    @test acceptance
    @printf("%17.12f | %17.12f %s\n", analytical, numerical, acceptance ? "✔" : "✗")
  #   end
  #   end
  # end
  # end
  # end
end

println("""```
""")

