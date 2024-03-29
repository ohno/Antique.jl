HO3D = HarmonicOscillator3D(k=1.0, m=1.0, ℏ=1.0)

# ∫|Rₙₗ(r)|²r²dr = δₙ₁ₙ₂δₗ₁ₗ₂


println(raw"""
#### Normalization of $R_{nl}(r)$

```math
\int |R_{nl}(r)|^2 r^2 \mathrm{d}r = 1
```
```""")

@testset "∫|Rₙₗ(r)|²r²dr = δₙ₁ₙ₂δₗ₁ₗ₂" begin
  println(" n |  l |        analytical |         numerical ")
  println("-- | -- | ----------------- | ----------------- ")
  for n in 0:5
  for l in 0:n
    analytical = 1
    numerical  = quadgk(r -> r^2 * R(HO3D,r,n=n,l=l)^2, 0, Inf, maxevals=10^3)[1]
    acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-5) : isapprox(analytical, numerical, rtol=1e-5)
    @test acceptance
    @printf("%2d | %2d | %17.12f | %17.12f %s\n", n, l, analytical, numerical, acceptance ? "✔" : "✗")
  end
  end
end

println("""```
""")


# <ψₙ₁ₗ₁ₘ₁|ψₙ₂ₗ₂ₘ₂> = δₙ₁ₙ₂δₗ₁ₗ₂δₘ₁ₘ₂


println(raw"""
#### Normalization & Orthogonality of $\psi_n(r,\theta,\varphi)$

```math
\int \psi_i^\ast(r,\theta,\varphi) \psi_j(r,\theta,\varphi) r^2 \sin(\theta) \mathrm{d}r \mathrm{d}\theta \mathrm{d}\varphi = \delta_{ij}
```
```""")


@testset "<ψₙ₁ₗ₁ₘ₁|ψₙ₂ₗ₂ₘ₂> = δₙ₁ₙ₂δₗ₁ₗ₂δₘ₁ₘ₂" begin
  println("n₁ | n₂ | l₁ | l₂ | m₁ | m₂ |        analytical |         numerical ")
  println("-- | -- | -- | -- | -- | -- | ----------------- | ----------------- ")
  for n1 in 0:1
  for n2 in 0:1
  for l1 in 0:n1
  for l2 in 0:n2
  for m1 in -l1:l1
  for m2 in -l2:l2
    analytical = (n1 == n2 ? 1 : 0) * (l1 == l2 ? 1 : 0) * (m1 == m2 ? 1 : 0)
    numerical = real(
      quadgk(phi ->
      quadgk(theta ->
      quadgk(r ->
        r^2 * sin(theta) * conj(ψ(HO3D,r,theta,phi,n=n1,l=l1,m=m1)) * ψ(HO3D,r,theta,phi,n=n2,l=l2,m=m2)
      , 0, Inf, maxevals=50)[1]
      , 0, π, maxevals=4)[1]
      , 0, 2π, maxevals=8)[1]
    )
    acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-2) : isapprox(analytical, numerical, rtol=1e-2)
    @test acceptance
    @printf("%2d | %2d | %2d | %2d | %2d | %2d | %17.12f | %17.12f %s\n", n1, n2, l1, l2, m1, m2, analytical, numerical, acceptance ? "✔" : "✗")
  end
  end
  end
  end
  end
  end
end

println("""```
""")
