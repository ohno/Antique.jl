RR = RigidRotor(m₁=1.0, m₂=1.0, R=1.0, ℏ=1.0)


# Pₙᵐ(x) = √(1-x²)ᵐ dᵐ/dxᵐ Pₙ(x); Pₙ(x) = 1/(2ⁿn!) dⁿ/dxⁿ (x²-1)ⁿ


println(raw"""
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

@testset "Pₙᵐ(x) = √(1-x²)ᵐ dᵐ/dxᵐ Pₙ(x); Pₙ(x) = 1/(2ⁿn!) dⁿ/dxⁿ (x²-1)ⁿ" begin
  for n in 0:4
  for m in 0:n
      # Rodrigues' formula
      @variables x
      Dn = n==0 ? x->x : Differential(x)^n         # dⁿ/dxⁿ
      Dm = m==0 ? x->x : Differential(x)^m         # dᵐ/dxᵐ
      a = 1 // (2^n * factorial(n))                # left
      b = (x^2 - 1)^n                              # right
      c = (1 - x^2)^(m//2) * Dm(a * Dn(b))         # Rodrigues' formula
      d = expand_derivatives(c)                    # expand dⁿ/dxⁿ and dᵐ/dxᵐ
      e = simplify(d, expand=true)                 # simplify
      f = simplify(P(RR, x, n=n, m=m), expand=true) # closed-form
      # latexify
      eq1 = latexify(e, env=:raw)
      eq2 = latexify(f, env=:raw)
      # judge
      acceptance = isequal(e, f)
      println("``n=$n, m=$m:`` ", acceptance ? "✔" : "✗")
      # show LaTeX
      println("""```math
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


println(raw"""
#### Normalization & Orthogonality of $P_n^m(x)$

```math
\int_{-1}^{1} P_i^m(x) P_j^m(x) \mathrm{d}x = \frac{2(j+m)!}{(2j+1)(j-m)!} \delta_{ij}
```
```""")

@testset "∫Pᵢᵐ(x)Pⱼᵐ(x)dx = 2(j+m)!/(2j+1)(j-m)! δᵢⱼ" begin
  println(" m |  i |  j |        analytical |         numerical ")
  println("-- | -- | -- | ----------------- | ----------------- ")
  for m in 0:5
  for i in m:9
  for j in m:9
    analytical = 2*factorial(j+m)/(2*j+1)/factorial(j-m)*(i == j ? 1 : 0)
    numerical  = quadgk(x -> P(RR, x, n=i, m=m) * P(RR, x, n=j, m=m), -1, 1, maxevals=10^3)[1]
    acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-5) : isapprox(analytical, numerical, rtol=1e-5)
    @test acceptance
    @printf("%2d | %2d | %2d | %17.12f | %17.12f %s\n", m, i, j, analytical, numerical, acceptance ? "✔" : "✗")
  end
  end
  end
end

println("""```
""")


# ∫Yₗ₁ₘ₁(θ,φ)Yₗ₂ₘ₂(θ,φ)sinθdθdφ = δₗ₁ₗ₂δₘ₁ₘ₂


println(raw"""
#### Normalization & Orthogonality of $Y_{lm}(\theta,\varphi)$

```math
\int_0^{2\pi}
\int_0^\pi
Y_{lm}(\theta,\varphi)^* Y_{l'm'}(\theta,\varphi) \sin(\theta)
~\mathrm{d}\theta \mathrm{d}\varphi
= \delta_{ll'} \delta_{mm'}
```
```""")

@testset "∫Yₗ₁ₘ₁(θ,φ)Yₗ₂ₘ₂(θ,φ)sinθdθdφ = δₗ₁ₗ₂δₘ₁ₘ₂" begin
  println("l₁ | l₂ | m₁ | m₂ |        analytical |         numerical ")
  println("-- | -- | -- | -- | ----------------- | ----------------- ")
  for l1 in 0:2
  for l2 in 0:2
  for m₁ in -l1:l1
  for m₂ in -l2:l2
    analytical = (l1 == l2 ? 1 : 0) * (m₁ == m₂ ? 1 : 0)
    numerical  = real(
      quadgk(φ ->
      quadgk(θ ->
        conj(Y(RR,θ,φ,l=l1,m=m₁)) * Y(RR,θ,φ,l=l2,m=m₂) * sin(θ)
      , 0, π, maxevals=50)[1]
      , 0, 2π, maxevals=100)[1]
    )
    acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-5) : isapprox(analytical, numerical, rtol=1e-5)
    @test acceptance
    @printf("%2d | %2d | %2d | %2d | %17.12f | %17.12f %s\n", l1, l2, m₁, m₂, analytical, numerical, acceptance ? "✔" : "✗")
  end
  end
  end
  end
end

println("""```
""")


# <ψₙ|H|ψₙ>  = ∫ψₙ*Tψₙdx = Eₙ

println(raw"""
#### Eigenvalues

```math
  \begin{aligned}
    E_n
    &=      \int_{0}^{2\pi}~\mathrm{d}\varphi \int_{0}^{\pi}~\mathrm{d}\theta  ~ \psi^\ast_{lm}(\theta,\varphi) \hat{H} \psi_{lm}(\theta,\varphi)  \\
    &=      \int_{0}^{2\pi}~\mathrm{d}\varphi \int_{0}^{\pi}~\mathrm{d}\theta  ~ \psi^\ast_{lm}(\theta,\varphi) \left[ \hat{V} + \hat{T} \right] \psi_{lm}(\theta,\varphi) \\
    &=      \int_{0}^{2\pi}~\mathrm{d}\varphi \int_{0}^{\pi}~\mathrm{d}\theta  ~ \psi^\ast_{lm}(\theta,\varphi) \left[ 0 + \frac{l(l+1)}{2I} \right] \psi_{lm}(\theta,\varphi) \\
    &=      \int_{0}^{2\pi}~\mathrm{d}\varphi \int_{0}^{\pi}~\mathrm{d}\theta  ~ \frac{l(l+1)}{2I} |\psi_{lm}(\theta,\varphi)|^2
  \end{aligned}
```
where $l$ is angular momentum quantum number and $I$ is the moment of intertia.
```""")

@testset "<ψₙ|H|ψₙ> = ∫ψₙ*Hψₙdx = Eₙ" begin
  function ψHψ(RR;l)
      μ = (1/RR.m₁ + 1/RR.m₂)^(-1)
      I = μ * RR.R^2
      YY = real(quadgk(φ -> quadgk(θ -> conj(Y(RR,θ,φ,l=l,m=0)) * Y(RR,θ,φ,l=l,m=0) * sin(θ), 0, π, maxevals=10)[1], 0, 2π, maxevals=10)[1])
      T = l*(l+1)/2/I * YY
      return T 
  end
  
  ℏ=1
  println(" m₁ |  m₂ |  R  |  l  |       analytical  |         numerical ")
  println("--- | --- | --- | --- | ----------------- | ----------------- ")
  for m₁ in [0.5, 2.0]
  for m₂ in [0.5, 2.0]
  for R in [0.5, 2.0]
  for l in [1, 2, 3]
      RR = RigidRotor(m₁=m₁, m₂=m₂, R=R, ℏ=ℏ)
      analytical = E(RR,l=l)
      numerical = ψHψ(RR,l=l)
      acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-4) : isapprox(analytical, numerical, rtol=1e-4)
      @printf("%.1f | %.1f | %.1f | %.1f | %17.12f | %17.12f %s\n", m₁, m₂, R, l,analytical, numerical, acceptance ? "✔" : "✗") 
  end
  end
  end
  end
end




println("""```
""")
