io = open("./result/SphericalOscillator.md", "w")
SO = SphericalOscillator(k=1.0, μ=1.0, ℏ=1.0)


# Pₙᵐ(x) = √(1-x²)ᵐ dᵐ/dxᵐ Pₙ(x); Pₙ(x) = 1/(2ⁿn!) dⁿ/dxⁿ (x²-1)ⁿ


println(io, raw"""
#### Associated Legendre Polynomials $P_n^m(x)$

```math
  \begin{aligned}
    P_n^m(x)
    &= \left( 1-x^2 \right)^{m/2} \frac{\mathrm{d}^m}{\mathrm{d}x^m} P_n(x) \\
    &= \left( 1-x^2 \right)^{m/2} \frac{\mathrm{d}^m}{\mathrm{d}x^m} \frac{1}{2^n n!} \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left[ \left( x^2-1 \right)^n \right] \\
    &= \frac{1}{2^n} (1-x^2)^{m/2} \sum_{j=0}^{\left\lfloor\frac{n-m}{2}\right\rfloor} (-1)^j \frac{(2n-2j)!}{j! (n-j)! (n-2j-m)!} x^{(n-2j-m)}.
  \end{aligned}
```
""")

@testset "SO: Pₙᵐ(x) = √(1-x²)ᵐ dᵐ/dxᵐ Pₙ(x); Pₙ(x) = 1/(2ⁿn!) dⁿ/dxⁿ (x²-1)ⁿ" begin
  @variables x
  for n in 0:4
  for m in 0:n
    # Rodrigues' formula
    Dn = n==0 ? x->x : Differential(x)^n                  # dⁿ/dxⁿ
    Dm = m==0 ? x->x : Differential(x)^m                  # dᵐ/dxᵐ
    a = 1 // (2^n * factorial(n))                         # left
    b = (x^2 - 1)^n                                       # right
    c = (1 - x^2)^(m//2) * Dm(a * Dn(b))                  # Rodrigues' formula
    d = expand_derivatives(c)                             # expand dⁿ/dxⁿ and dᵐ/dxᵐ
    e = simplify(d, expand=true)                          # simplify
    f = simplify(Antique.P(SO, x, n=n, m=m), expand=true) # closed-form
    # latexify
    eq1 = latexify(e, env=:raw)
    eq2 = latexify(f, env=:raw)
    # judge
    acceptance = isequal(e, f)
    println(io, "``n=$n, m=$m:`` ", acceptance ? "✔" : "✗")
    # show LaTeX
    println(io, """```math
    \\begin{aligned}
      P_{$n}^{$m}(x)
       = $(latexify(c, env=:raw))
      &= $(eq1) \\\\
      &= $(eq2)
    \\end{aligned}
    ```
    """)
    # result
    @test acceptance
  end
  end
end


# ∫Pᵢᵐ(x)Pⱼᵐ(x)dx = 2(j+m)!/(2j+1)(j-m)! δᵢⱼ


println(io, raw"""
#### Normalization & Orthogonality of $P_n^m(x)$

```math
\int_{-1}^{1} P_i^m(x) P_j^m(x) \mathrm{d}x = \frac{2(j+m)!}{(2j+1)(j-m)!} \delta_{ij}
```

```""")

@testset "SO: ∫Pᵢᵐ(x)Pⱼᵐ(x)dx = 2(j+m)!/(2j+1)(j-m)! δᵢⱼ" begin
  println(io, " m |  i |  j |     analytical |      numerical ")
  println(io, "-- | -- | -- | -------------- | -------------- ")
  for m in 0:5
  for i in m:9
  for j in m:9
    analytical = 2*factorial(j+m)/(2*j+1)/factorial(j-m)*(i == j ? 1 : 0)
    numerical  = quadgk(x -> Antique.P(SO, x, n=i, m=m) * Antique.P(SO, x, n=j, m=m), -1, 1, maxevals=10^3)[1]
    acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-5) : isapprox(analytical, numerical, rtol=1e-5)
    @test acceptance
    @printf(io, "%2d | %2d | %2d | %14.9f | %14.9f %s\n", m, i, j, analytical, numerical, acceptance ? "✔" : "✗")
  end
  end
  end
end

println(io, """```\n""")


# ∫Yₗ₁ₘ₁(θ,φ)Yₗ₂ₘ₂(θ,φ)sinθdθdφ = δₗ₁ₗ₂δₘ₁ₘ₂


println(io, raw"""
#### Normalization & Orthogonality of $Y_{lm}(\theta,\varphi)$

```math
\int_0^{2\pi}
\int_0^\pi
Y_{lm}(\theta,\varphi)^* Y_{l'm'}(\theta,\varphi) \sin(\theta)
~\mathrm{d}\theta \mathrm{d}\varphi
= \delta_{ll'} \delta_{mm'}
```

```""")

@testset "SO: ∫Yₗ₁ₘ₁(θ,φ)Yₗ₂ₘ₂(θ,φ)sinθdθdφ = δₗ₁ₗ₂δₘ₁ₘ₂" begin
  println(io, "l₁ | l₂ | m₁ | m₂ |     analytical |      numerical ")
  println(io, "-- | -- | -- | -- | -------------- | -------------- ")
  for l1 in 0:2
  for l2 in 0:2
  for m1 in -l1:l1
  for m2 in -l2:l2
    analytical = (l1 == l2 ? 1 : 0) * (m1 == m2 ? 1 : 0)
    numerical  = real(
      quadgk(φ ->
      quadgk(θ ->
        conj(Antique.Y(SO,θ,φ,l=l1,m=m1)) * Antique.Y(SO,θ,φ,l=l2,m=m2) * sin(θ)
      , 0, π, maxevals=50)[1]
      , 0, 2π, maxevals=100)[1]
    )
    acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-5) : isapprox(analytical, numerical, rtol=1e-5)
    @test acceptance
    @printf(io, "%2d | %2d | %2d | %2d | %14.9f | %14.9f %s\n", l1, l2, m1, m2, analytical, numerical, acceptance ? "✔" : "✗")
  end
  end
  end
  end
end

println(io, """```\n""")


# Lₙ⁽ᵅ⁾(x) = x⁻ᵅeˣ/n! dⁿ/dxⁿ xⁿ⁺ᵅe⁻ˣ


println(io, raw"""
#### Generalized Laguerre Polynomials $L_n^{(\alpha)}(x)$

```math
  \begin{aligned}
    L_n^{(\alpha)}(x)
    &= \frac{x^{-\alpha}e^x}{n!} \frac{d^n}{dx^n}\left(x^{n+\alpha}e^{-x}\right) \\
    &= \sum_{k=0}^n(-1)^k \frac{\Gamma(\alpha+n+1)}{\Gamma(\alpha+k+1)\Gamma(n-k+1)} \frac{x^k}{k !}.
  \end{aligned}
```
""")

@testset "SO: Lₙ⁽ᵅ⁾(x) = x⁻ᵅeˣ/n! dⁿ/dxⁿ xⁿ⁺ᵅe⁻ˣ" begin
  @variables x
  for n in 0:4
  for α in 0:n
    # Rodrigues' formula
    D = n==0 ? x->x : Differential(x)^n                   # dⁿ/dxⁿ
    a = exp(x) * x^(-α) / factorial(n)                    # left
    b = exp(-x) * x^(n+α)                                 # right
    c = a * D(b)                                          # Rodrigues' formula
    d = expand_derivatives(c)                             # expand dⁿ/dxⁿ
    e = simplify(d, expand=true)                          # simplify
    f = simplify(Antique.L(SO, x, n=n, α=α), expand=true) # closed-form
    # latexify
    eq1 = latexify(e, env=:raw)
    eq2 = latexify(f, env=:raw)
    # judge
    acceptance = isequal(e, f)
    println(io, "``n=$n, α=$α:`` ", acceptance ? "✔" : "✗")
    # show LaTeX
    println(io, """```math
    \\begin{aligned}
      L_{$n}^{($α)}(x)
       = $(latexify(c, env=:raw))
      &= $(eq1) \\\\
      &= $(eq2)
    \\end{aligned}
    ```
    """)
    # result
    @test acceptance
  end
  end
end


# ∫Lᵢ⁽ᵅ⁾(x)Lⱼ⁽ᵅ⁾(x)xᵅexp(-x)dx = Γ(i+α+1)/i! δᵢⱼ


println(io, raw"""
#### Normalization & Orthogonality of $L_n^{(\alpha)}(x)$

```math
\int_0^\infty L_i^{(\alpha)}(x) L_j^{(\alpha)}(x) x^\alpha \mathrm{e}^{-x} \mathrm{d}x = \frac{\Gamma(n+\alpha+1)}{n!} \delta_{ij}
```

```""")

@testset "SO: ∫Lᵢ⁽ᵅ⁾(x)Lⱼ⁽ᵅ⁾(x)xᵅexp(-x)dx = Γ(i+α+1)/i! δᵢⱼ" begin
  println(io, "   α |  i |  j |     analytical |      numerical ")
  println(io, "---- | -- | -- | -------------- | -------------- ")
  for α in [0.01, 0.1, 1.0]
  for i in 0:9
  for j in 0:9
    analytical = gamma(i+α+1)/factorial(i)*(i == j ? 1 : 0)
    numerical  = quadgk(x -> real(Antique.L(SO, x, n=i, α=α)) * real(Antique.L(SO, x, n=j, α=α)) * x^α * exp(-x), 0, Inf, maxevals=10^3)[1]
    acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-5) : isapprox(analytical, numerical, rtol=1e-5)
    @test acceptance
    @printf(io, "%4.2f | %2d | %2d | %14.9f | %14.9f %s\n", α, i, j, analytical, numerical, acceptance ? "✔" : "✗")
  end
  end
  end
end

println(io, """```\n""")


# ∫|Rₙₗ(r)|²r²dr = δₙ₁ₙ₂δₗ₁ₗ₂


println(io, raw"""
#### Normalization of $R_{nl}(r)$

```math
\int |R_{nl}(r)|^2 r^2 \mathrm{d}r = 1
```

```""")

@testset "SO: ∫|Rₙₗ(r)|²r²dr = δₙ₁ₙ₂δₗ₁ₗ₂" begin
  println(io, " n |  l |     analytical |      numerical ")
  println(io, "-- | -- | -------------- | -------------- ")
  for n in 0:5
  for l in 0:n
    analytical = 1
    numerical  = quadgk(r ->Antique.R(SO,r,n=n,l=l)^2 * r^2, 0, Inf, maxevals=10^3)[1]
    acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-5) : isapprox(analytical, numerical, rtol=1e-5)
    @test acceptance
    @printf(io, "%2d | %2d | %14.9f | %14.9f %s\n", n, l, analytical, numerical, acceptance ? "✔" : "✗")
  end
  end
end

println(io, """```\n""")


# 2 * <ψₙ|V|ψₙ> = Eₙ


println(io, raw"""
#### Virial Theorem

The virial theorem $\langle T \rangle = \langle V \rangle$ and the definition of Hamiltonian $\langle H \rangle = \langle T \rangle + \langle V \rangle$ derive $\langle H \rangle = 2 \langle V \rangle$.

```math
2 \langle V \rangle = 2 \times \int \psi_i^\ast(r,\theta,\varphi) V(r) \psi_j(r,\theta,\varphi) r^2 \sin(\theta) \mathrm{d}r \mathrm{d}\theta \mathrm{d}\varphi = 2 \times \int V(r) |R_{nl}(r)|^2 r^2 \mathrm{d}r = E_n
```

```""")

@testset "SO: 2 * <ψₙₗₘ|V|ψₙₗₘ> = 2∫ψ*ₙₗₘ(r,θ,φ)V(r)ψₙₗₘ(r,θ,φ)r²sin(θ)drdθdφ = Eₙ" begin
  println(io, " n |  l |  m |     analytical |      numerical ")
  println(io, "-- | -- | -- | -------------- | -------------- ")
  for n in 0:3
  for l in 0:n
  for m in -l:l
    analytical = E(SO, n=n, l=l)
    # numerical = real(
    #   quadgk(phi ->
    #   quadgk(theta ->
    #   quadgk(r ->
    #     2 * V(SO,r) * r^2 * sin(theta) * conj(ψ(SO,r,theta,phi,n=n,l=l,m=m)) * ψ(SO,r,theta,phi,n=n,l=l,m=m)
    #   , 0, Inf, maxevals=50)[1]
    #   , 0, π, maxevals=4)[1]
    #   , 0, 2π, maxevals=8)[1]
    # )
    numerical = real(first(hcubature(r -> 2 * V(SO,r[1]) * r[1]^2 * sin(r[2]) * conj(ψ(SO,r[1],r[2],r[3],n=n,l=l,m=m)) * ψ(SO,r[1],r[2],r[3],n=n,l=l,m=m), [0,0,0], [100,π,2π], maxevals=2000)))
    acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-1) : isapprox(analytical, numerical, rtol=1e-1)
    @test acceptance
    @printf(io, "%2d | %2d | %2d | %14.9f | %14.9f %s\n", n, l, m, analytical, numerical, acceptance ? "✔" : "✗")
  end
  end
  end
end

println(io, """```

```""")

@testset "SO: 2 * <ψₙₗₘ|V|ψₙₗₘ> = 2∫V(r)|Rₙₗ(r)|²r²dr = Eₙ" begin
  println(io, " n |  l |     analytical |      numerical ")
  println(io, "-- | -- | -------------- | -------------- ")
  for n in 0:5
    for l in 0:n
      analytical = E(SO, n=n, l=l)
      numerical  = 2 * quadgk(r -> V(SO,r) *Antique.R(SO,r,n=n,l=l)^2 * r^2, 0, Inf, maxevals=10^3)[1]
      acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-5) : isapprox(analytical, numerical, rtol=1e-5)
      @test acceptance
      @printf(io, "%2d | %2d | %14.9f | %14.9f %s\n", n, l, analytical, numerical, acceptance ? "✔" : "✗")
    end
  end
end

println(io, """```\n""")


# <ψₙ₁ₗ₁ₘ₁|ψₙ₂ₗ₂ₘ₂> = δₙ₁ₙ₂δₗ₁ₗ₂δₘ₁ₘ₂


println(io, raw"""
#### Normalization & Orthogonality of $\psi_n(r,\theta,\varphi)$

```math
\int \psi_i^\ast(r,\theta,\varphi) \psi_j(r,\theta,\varphi) r^2 \sin(\theta) \mathrm{d}r \mathrm{d}\theta \mathrm{d}\varphi = \delta_{ij}
```

```""")

@testset "SO: <ψₙ₁ₗ₁ₘ₁|ψₙ₂ₗ₂ₘ₂> = δₙ₁ₙ₂δₗ₁ₗ₂δₘ₁ₘ₂" begin
  println(io, "n₁ | n₂ | l₁ | l₂ | m₁ | m₂ |     analytical |      numerical ")
  println(io, "-- | -- | -- | -- | -- | -- | -------------- | -------------- ")
  for n1 in 0:2
  for n2 in 0:2
  for l1 in 0:n1
  for l2 in 0:n2
  for m1 in -l1:l1
  for m2 in -l2:l2
    analytical = (n1 == n2 ? 1 : 0) * (l1 == l2 ? 1 : 0) * (m1 == m2 ? 1 : 0)
    # numerical = real(
    #   quadgk(phi ->
    #   quadgk(theta ->
    #   quadgk(r ->
    #     r^2 * sin(theta) * conj(ψ(SO,r,theta,phi,n=n1,l=l1,m=m1)) * ψ(SO,r,theta,phi,n=n2,l=l2,m=m2)
    #   , 0, Inf, maxevals=50)[1]
    #   , 0, π, maxevals=4)[1]
    #   , 0, 2π, maxevals=8)[1]
    # )
    numerical = real(first(hcubature(r -> r[1]^2 * sin(r[2]) * conj(ψ(SO,r[1],r[2],r[3],n=n1,l=l1,m=m1)) * ψ(SO,r[1],r[2],r[3],n=n2,l=l2,m=m2), [0,0,0], [100,π,2π], maxevals=2000)))
    acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-1) : isapprox(analytical, numerical, rtol=1e-1)
    @test acceptance
    @printf(io, "%2d | %2d | %2d | %2d | %2d | %2d | %14.9f | %14.9f %s\n", n1, n2, l1, l2, m1, m2, analytical, numerical, acceptance ? "✔" : "✗")
  end
  end
  end
  end
  end
  end
end

println(io, """```\n""")


close(io)