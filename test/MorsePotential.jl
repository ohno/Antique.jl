MP = MorsePotential()


# Lₙ⁽ᵅ⁾(x) = x⁻ᵅeˣ/n! dⁿ/dxⁿ xⁿ⁺ᵅe⁻ˣ


println(raw"""
#### Generalized Laguerre Polynomials $L_n^{(\alpha)}(x)$

```math
  \begin{aligned}
    L_n^{(\alpha)}(x)
    &= \frac{x^{-\alpha}e^x}{n!} \frac{d^n}{dx^n}\left(x^{n+\alpha}e^{-x}\right) \\
    &= \sum_{k=0}^n(-1)^k \frac{\Gamma(\alpha+n+1)}{\Gamma(\alpha+k+1)\Gamma(n-k+1)} \frac{x^k}{k !}.
  \end{aligned}
```
""")

@testset "Lₙ⁽ᵅ⁾(x) = x⁻ᵅeˣ/n! dⁿ/dxⁿ xⁿ⁺ᵅe⁻ˣ" begin
  for n in 0:4
  for α in 0:n
    # Rodriguesの公式の展開
    @variables x
    D = n==0 ? x->x : Differential(x)^n           # dⁿ/dxⁿ
    a = exp(x) * x^(-α) / factorial(n)            # left
    b = exp(-x) * x^(n+α)                         # right
    c = a * D(b)                                  # Rodrigues' formula
    d = expand_derivatives(c)                     # expand dⁿ/dxⁿ
    e = simplify(d, expand=true)                  # simplify
    f = simplify(L(MP, x, n=n, α=α), expand=true) # closed-form
    # latexify
    eq1 = latexify(e, env=:raw)
    eq2 = latexify(f, env=:raw)
    # judge
    acceptance = isequal(e, f)
    println("``n=$n, α=$α:`` ", acceptance ? "✔" : "✗")
    # show LaTeX
    println("""```math
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


println(raw"""
#### Normalization & Orthogonality of $L_n^{(\alpha)}(x)$

```math
\int_0^\infty L_i^{(\alpha)}(x) L_j^{(\alpha)}(x) x^\alpha \mathrm{e}^{-x} \mathrm{d}x = \frac{\Gamma(n+\alpha+1)}{n!} \delta_{ij}
```

```""")

@testset "∫Lᵢ⁽ᵅ⁾(x)Lⱼ⁽ᵅ⁾(x)xᵅexp(-x)dx = Γ(i+α+1)/i! δᵢⱼ" begin
  println("   α |  i |  j |        analytical |         numerical ")
  println("---- | -- | -- | ----------------- | ----------------- ")
  for α in [0.01,0.05,0.1,0.5,1.0]
  for i in 0:9
  for j in 0:9
    analytical = gamma(i+α+1)/factorial(i)*(i == j ? 1 : 0)
    numerical  = quadgk(x -> real(L(MP, x, n=i, α=α)) * real(L(MP, x, n=j, α=α)) * x^α * exp(-x), 0, Inf, maxevals=10^3)[1]
    acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-5) : isapprox(analytical, numerical, rtol=1e-5)
    @test acceptance
    @printf("%4.2f | %2d | %2d | %17.12f | %17.12f %s\n", α, i, j, analytical, numerical, acceptance ? "✔" : "✗")
  end
  end
  end
end

println("""```
""")


# <ψᵢ|ψⱼ> = δᵢⱼ


println(raw"""
#### Normalization & Orthogonality of $\psi_n(r)$

```math
\int_0^\infty \psi_i^\ast(r) \psi_j(r) \mathrm{d}r = \delta_{ij}
```

```""")

@testset "<ψᵢ|ψⱼ> = δᵢⱼ" begin
  println(" i |  j |        analytical |         numerical ")
  println("-- | -- | ----------------- | ----------------- ")
  for i in 0:9
  for j in 0:9
    analytical = (i == j ? 1 : 0)
    numerical  = quadgk(x -> conj(ψ(MP, x, n=i)) * ψ(MP, x, n=j), 0, Inf, maxevals=10^3)[1]
    acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-5) : isapprox(analytical, numerical, rtol=1e-5)
    @test acceptance
    @printf("%2d | %2d | %17.12f | %17.12f %s\n", i, j, analytical, numerical, acceptance ? "✔" : "✗")
  end
  end
end

println("""```
""")


# <ψₙ|H|ψₙ> = ∫ψₙ*Hψₙdx = Eₙ


println(raw"""
#### Eigenvalues

```math
  \begin{aligned}
    E_n
    &=      \int \psi^\ast_n(r) \hat{H} \psi_n(r) \mathrm{d}r \\
    &=      \int \psi^\ast_n(r) \left[ \hat{V} + \hat{T} \right] \psi(r) \mathrm{d}r \\
    &=      \int \psi^\ast_n(r) \left[ V(r) - \frac{\hbar^2}{2m} \frac{\mathrm{d}^{2}}{\mathrm{d} r^{2}} \right] \psi(r) \mathrm{d}r \\
    &\simeq \int \psi^\ast_n(r) \left[ V(r)\psi(r) -\frac{\hbar^2}{2m} \frac{\psi(r+\Delta r) - 2\psi(r) + \psi(r-\Delta r)}{\Delta r^{2}} \right] \mathrm{d}r.
  \end{aligned}
```

Where, the difference formula for the 2nd-order derivative:

```math
\begin{aligned}
  % 2\psi(r)
  % + \frac{\mathrm{d}^{2} \psi(r)}{\mathrm{d} r^{2}} \Delta r^{2}
  % + O\left(\Delta r^{4}\right)
  % &=
  % \psi(r+\Delta r)
  % + \psi(r-\Delta r)
  % \\
  % \frac{\mathrm{d}^{2} \psi(r)}{\mathrm{d} r^{2}} \Delta r^{2}
  % &=
  % \psi(r+\Delta r)
  % - 2\psi(r)
  % + \psi(r-\Delta r)
  % - O\left(\Delta r^{4}\right)
  % \\
  % \frac{\mathrm{d}^{2} \psi(r)}{\mathrm{d} r^{2}}
  % &=
  % \frac{\psi(r+\Delta r) - 2\psi(r) + \psi(r-\Delta r)}{\Delta r^{2}}
  % - \frac{O\left(\Delta r^{4}\right)}{\Delta r^{2}}
  % \\
  \frac{\mathrm{d}^{2} \psi(r)}{\mathrm{d} r^{2}}
  &=
  \frac{\psi(r+\Delta r) - 2\psi(r) + \psi(r-\Delta r)}{\Delta r^{2}}
  + O\left(\Delta r^{2}\right)
\end{aligned}
```

are given by the sum of 2 Taylor series:

```math
\begin{aligned}
\psi(r+\Delta r)
&= \psi(r)
+ \frac{\mathrm{d} \psi(r)}{\mathrm{d} r} \Delta r
+ \frac{1}{2!} \frac{\mathrm{d}^{2} \psi(r)}{\mathrm{d} r^{2}} \Delta r^{2}
+ \frac{1}{3!} \frac{\mathrm{d}^{3} \psi(r)}{\mathrm{d} r^{3}} \Delta r^{3}
+ O\left(\Delta r^{4}\right),
\\
\psi(r-\Delta r)
&= \psi(r)
- \frac{\mathrm{d} \psi(r)}{\mathrm{d} r} \Delta r
+ \frac{1}{2!} \frac{\mathrm{d}^{2} \psi(r)}{\mathrm{d} r^{2}} \Delta r^{2}
- \frac{1}{3!} \frac{\mathrm{d}^{3} \psi(r)}{\mathrm{d} r^{3}} \Delta r^{3}
+ O\left(\Delta r^{4}\right).
\end{aligned}
```

```""")

@testset "<ψₙ|H|ψₙ> = ∫ψₙ*Hψₙdx = Eₙ" begin
  ψHψ(MP, r; n=0, Δr=0.005) = V(MP,r)*ψ(MP,r,n=n)^2 - MP.ℏ^2/(2*MP.μ)*conj(ψ(MP,r,n=n))*(ψ(MP,r+Δr,n=n)-2*ψ(MP,r,n=n)+ψ(MP,r-Δr,n=n))/Δr^2
  println("  k |  n |        analytical |         numerical ")
  println("--- | -- | ----------------- | ----------------- ")
  for k in [0.1,0.2,0.3,MP.k]
  for n in 0:9
    MP = MorsePotential(k=k)
    analytical = E(MP, n=n)
    numerical  = quadgk(r -> ψHψ(MP, r, n=n, Δr=0.0001), 0.0001, Inf, maxevals=10^4)[1]
    acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-5) : isapprox(analytical, numerical, rtol=1e-5)
    @test acceptance
    @printf("%.1f | %2d | %17.12f | %17.12f %s\n", k, n, analytical, numerical, acceptance ? "✔" : "✗")
  end
  end
end

println("""```
""")


# 0 < Eₙ₊₁ - Eₙ for 0 ≤ n ≤ nₘₐₓ


println(raw"""
#### Recurrence Relation between $E_{n+1}$ and $E_n$

```math
\begin{equation}
\left\{ \,
  \begin{aligned}
    0 < \Delta E && 0 \leq n \leq n_\mathrm{max} \\
    \Delta E < 0 && \mathrm{otherwise}
  \end{aligned}
\right.
\end{equation}
```

```math
\Delta E =  E_{n+1} - E_n
```

```math
n_\mathrm{max} = \left\lfloor\frac{2 D_{\mathrm{e}}-h \nu_0}{h \nu_0}\right\rfloor
```

```""")

@testset "0 < Eₙ₊₁ - Eₙ for 0 ≤ n ≤ nₘₐₓ" begin
  println(" n  Eₙ          ΔE")
  for n in 0:nₘₐₓ(MP)+5
    Eₙ₊₁ = E(MP, n=n+1, nocheck=true)
    Eₙ   = E(MP, n=n, nocheck=true)
    ΔE  = Eₙ₊₁ - Eₙ
    @printf("%2d  %+9.6f  %+9.6f  ", n, Eₙ, ΔE)
    if n ≤ nₘₐₓ(MP)
      print("0 < ΔE  ")
      acceptance = 0 < ΔE
    else
      print("ΔE < 0  ")
      acceptance = ΔE < 0
    end
    @test acceptance
    println(acceptance ? "✔" : "✗")
    if nₘₐₓ(MP) == n
      println("-----------------------------------  nₘₐₓ(MP) = ", nₘₐₓ(MP))
    end
  end
end

println("""```
""")

