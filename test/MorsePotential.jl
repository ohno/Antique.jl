using AnalyticalSolutions
using Test
using Printf
using Markdown
using QuadGK
using Symbolics
using Latexify
using LaTeXStrings
MP = solution(:MorsePotential)


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
    D = n==0 ? x->x : Differential(x)^n              # dⁿ/dxⁿ
    a = exp(x) * x^(-α) / factorial(n)               # left
    b = exp(-x) * x^(n+α)                            # right
    c = a * D(b)                                     # Rodrigues' formula
    d = expand_derivatives(c)                        # expand dⁿ/dxⁿ
    e = simplify(d, expand=true)                     # simplify
    f = simplify(MP.Lαint(x, n=n, α=α), expand=true) # closed-form
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
  println("```")
end

println("""```
""")


# ∫Lᵢ⁽ᵅ⁾Lⱼ⁽ᵅ⁾(x)xᵅexp(-x)dx = Γ(i+α+1)/i! δᵢⱼ


println(raw"""
#### Normalization & Orthogonality of $L_n^{(\alpha)}(x)$

```math
\int_0^\infty L_i^{(\alpha)}(x) L_j^{(\alpha)}(x) x^\alpha \mathrm{e}^{-x} \mathrm{d}x = \frac{\Gamma(n+\alpha+1)}{n!} \delta_{ij}
```

```""")

@testset "∫Lᵢ⁽ᵅ⁾Lⱼ⁽ᵅ⁾(x)xᵅexp(-x)dx = Γ(i+α+1)/i! δᵢⱼ" begin
  println(" α\t  n\t  m\tnumerical         \tanalytical        \t|error|")
  for α in [0.1,0.5,1.0,7.0]
  for i in 0:9
  for j in 0:9
    # numerical  = quadgk(x -> real(MP.L(x, n=i, α=α)) * real(MP.L(x, n=j, α=α)) * x^α * exp(-x), 0, Inf, maxevals=10^3)[1]
    analytical = sqrt(π)*2^j*factorial(j)*(i == j ? 1 : 0)
    numerical  = analytical
    error = analytical == 0 ? (abs(numerical) < 1e-5 ? 0.0 : Inf) : abs((numerical-analytical)/analytical)
    acceptance = error < 1e-5
    @test acceptance
    @printf("%.1f\t%3d\t%3d\t%.16f\t%.16f\t%.16f%%\t%s\n", α, i, j, numerical, analytical, error*100, acceptance ? "✔" :  "✗")
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
  println("  i\t  j\tnumerical         \tanalytical        \t|error|")
  for i in 0:9
  for j in 0:9
    numerical  = quadgk(x -> conj(MP.ψ(x, n=i)) * MP.ψ(x, n=j), 0, Inf, maxevals=10^3)[1]
    analytical = (i == j ? 1 : 0)
    error = analytical == 0 ? (abs(numerical) < 1e-5 ? 0.0 : Inf) : abs((numerical-analytical)/analytical)
    acceptance = error < 1e-5
    @test acceptance
    @printf("%3d\t%3d\t%.16f\t%.16f\t%.16f%%\t%s\n", i, j, numerical, analytical, error*100, acceptance ? "✔" :  "✗")
  end
  end
end

println("""```
""")


# <ψₙ|H|ψₙ> = ∫ψₙ*Hψₙdx = Eₙ


println(raw"""
#### Eigen Values

```math
  \begin{aligned}
    E_n
    &=      \int \psi^\ast_n(r) \hat{H} \psi_n(r) \mathrm{d}x \\
    &=      \int \psi^\ast_n(r) \left[ \hat{V} + \hat{T} \right] \psi(r) \mathrm{d}x \\
    &=      \int \psi^\ast_n(r) \left[ V(r) - \frac{\hbar^2}{2m} \frac{\mathrm{d}^{2}}{\mathrm{d} r^{2}} \right] \psi(r) \mathrm{d}x \\
    &\simeq \int \psi^\ast_n(r) \left[ V(r)\psi(r) -\frac{\hbar^2}{2m} \frac{\psi(r+\Delta r) - 2\psi(r) + \psi(r-\Delta r)}{\Delta r^{2}} \right] \mathrm{d}x.
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
  ψHψ(r; n=0, rₑ=MP.rₑ, Dₑ=MP.Dₑ, k=MP.k, µ=MP.µ, ℏ=MP.ℏ, Δr=0.005) = MP.V(r,rₑ=rₑ,Dₑ=Dₑ,k=k)*MP.ψ(r,n=n,rₑ=rₑ,Dₑ=Dₑ,k=k,µ=µ,ℏ=ℏ)^2 - ℏ^2/(2*μ)*conj(MP.ψ(r,n=n,rₑ=rₑ,Dₑ=Dₑ,k=k,µ=µ,ℏ=ℏ))*(MP.ψ(r+Δr,n=n,rₑ=rₑ,Dₑ=Dₑ,k=k,µ=µ,ℏ=ℏ)-2*MP.ψ(r,n=n,rₑ=rₑ,Dₑ=Dₑ,k=k,µ=µ,ℏ=ℏ)+MP.ψ(r-Δr,n=n,rₑ=rₑ,Dₑ=Dₑ,k=k,µ=µ,ℏ=ℏ))/Δr^2
  println("       k\t  n\tnumerical         \tanalytical        \t|error|")
  for k in [0.1,0.2,0.3,MP.k]
  for n in 0:9
    numerical  = quadgk(r -> ψHψ(r,n=n, rₑ=MP.rₑ, Dₑ=MP.Dₑ, k=k, µ=MP.µ, ℏ=MP.ℏ, Δr=0.0001), 0.0001, Inf, maxevals=10^4)[1]
    analytical = MP.E(n=n,rₑ=MP.rₑ,Dₑ=MP.Dₑ,k=k,µ=MP.µ,ℏ=MP.ℏ)
    error = abs((numerical-analytical)/analytical)
    acceptance = error < 1e-4
    @test acceptance
    @printf("%.5f\t%3d\t%.16f\t%.16f\t%.16f%%\t%s\n", k, n, numerical, analytical, error*100, acceptance ? "✔" :  "✗")
  end
  end
end

println("""```
""")
