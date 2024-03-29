PT = PoschlTeller(λ=10.0)

# <ψᵢ|ψⱼ> = δᵢⱼ


println(raw"""
#### Normalization & Orthogonality of $\psi_n(x)$

```math
\int \psi_i^\ast(x) \psi_j(x) \mathrm{d}x = \delta_{ij}
```

```""")

@testset "<ψᵢ|ψⱼ> = δᵢⱼ" begin
  println(" i |  j |        analytical |         numerical ")
  println("-- | -- | ----------------- | ----------------- ")
  for i in 0:9
  for j in 0:9
    analytical = (i == j ? 1 : 0)
    numerical  = quadgk(x -> conj(ψ(PT, x, n=i)) * ψ(PT, x, n=j), -Inf, Inf, maxevals=10^3)[1]
    acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-5) : isapprox(analytical, numerical, rtol=1e-5)
    @test acceptance
    @printf("%2d | %2d | %17.12f | %17.12f %s\n", i, j, analytical, numerical, acceptance ? "✔" : "✗")
  end
  end
end

println("""```
""")