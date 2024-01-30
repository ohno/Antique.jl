using Antique
using Test
using Printf
using Markdown
using QuadGK
using Symbolics
using Latexify
using LaTeXStrings
HO = HarmonicOscillator(k=1.0, m=1.0, ℏ=1.0)


# Hₙ(x) = (-1)ⁿ exp(x²) dⁿ/dxⁿ　exp(-x²) = ...


println(raw"""
#### Hermite Polynomials $H_n(x)$

```math
  \begin{aligned}
    H_{n}(x)
    &:= (-1)^n \mathrm{e}^{x^2} \frac{\mathrm{d}^n}{\mathrm{d}x^n} \mathrm{e}^{-x^2} \\
    &= n! \sum_{m=0}^{\lfloor n/2 \rfloor} \frac{(-1)^m}{m! (n-2m)!}(2 x)^{n-2m}.
  \end{aligned}
```
""")

@testset "Hₙ(x) = (-1)ⁿ exp(x²) dⁿ/dxⁿ exp(-x²) = ..." begin
  for n in 0:9
    # Rodrigues' formula
    @variables x
    D = n==0 ? x->x : Differential(x)^n     # dⁿ/dxⁿ
    a = (-1)^n * exp(x^2)                   # left
    b = exp(-x^2)                           # right
    c = a * D(b)                            # Rodrigues' formula
    d = expand_derivatives(c)               # expand dⁿ/dxⁿ
    e = simplify(d, expand=true)            # simplify
    f = simplify(H(HO, x, n=n), expand=true) # closed-form
    # latexify
    eq1 = latexify(e, env=:raw)
    eq2 = latexify(f, env=:raw)
    # judge
    acceptance = isequal(e, f)
    println("``n=$n:`` ", acceptance ? "✔" : "✗")
    # show LaTeX
    println("""```math
    \\begin{aligned}
      H_{$n}(x)
       = $(latexify(c, env=:raw))
      &= $(eq1) \\\\
      &= $(eq2)
    \\end{aligned}
    ```
    """)
    # result
    @test acceptance
  end
  println("```")
end

println("""```
""")


# ∫Hⱼ(x)Hᵢ(x)exp(-x²)dx = √π2ʲj!δᵢⱼ


println(raw"""
#### Normalization & Orthogonality of $H_n(x)$

```math
\int_{-\infty}^\infty H_j(x) H_i(x) \mathrm{e}^{-x^2} \mathrm{d}x = \sqrt{\pi} 2^j j! \delta_{ij}
```

```""")

@testset "∫Hⱼ(x)Hᵢ(x)exp(-x²)dx = √π2ʲj!δᵢⱼ" begin
  println(" i |  j |        analytical |         numerical ")
  println("-- | -- | ----------------- | ----------------- ")
  for i in 0:9
  for j in 0:9
    analytical = sqrt(π)*2^j*factorial(j)*(i == j ? 1 : 0)
    numerical  = quadgk(x -> H(HO, x, n=j) * H(HO, x, n=i)* exp(-x^2), -Inf, Inf, maxevals=10^3)[1]
    acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-5) : isapprox(analytical, numerical, rtol=1e-5)
    @test acceptance
    @printf("%2d | %2d | %17.12f | %17.12f %s\n", i, j, analytical, numerical, acceptance ? "✔" : "✗")
  end
  end
end

println("""```
""")


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
    numerical  = quadgk(x -> conj(ψ(HO, x, n=i)) * ψ(HO, x, n=j), -Inf, Inf, maxevals=10^3)[1]
    acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-5) : isapprox(analytical, numerical, rtol=1e-5)
    @test acceptance
    @printf("%2d | %2d | %17.12f | %17.12f %s\n", i, j, analytical, numerical, acceptance ? "✔" : "✗")
  end
  end
end

println("""```
""")


# 2 × <ψₙ|V|ψₙ> = Eₙ


println(raw"""
#### Virial Theorem

The virial theorem $\langle T \rangle = \langle V \rangle$ and the definition of Hamiltonian $\langle H \rangle = \langle T \rangle + \langle V \rangle$ derive $\langle H \rangle = 2 \langle V \rangle = 2 \langle T \rangle$.

```math
2 \int \psi_n^\ast(x) V(x) \psi_n(x) \mathrm{d}x = E_n
```

```""")

@testset "2 × <ψₙ|V|ψₙ> = Eₙ" begin
  println("  k |  n |        analytical |         numerical ")
  println("--- | -- | ----------------- | ----------------- ")
  for k in [0.1,0.5,1.0,5.0]
  for n in 0:9
    analytical = E(HO, n=n)
    numerical  = 2 * quadgk(x -> conj(ψ(HO, x, n=n)) * V(HO, x) * ψ(HO, x, n=n), -Inf, Inf, maxevals=10^3)[1]
    acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-5) : isapprox(analytical, numerical, rtol=1e-5)
    @test acceptance
    @printf("%.1f | %2d | %17.12f | %17.12f %s\n", k, n, analytical, numerical, acceptance ? "✔" : "✗")
  end
  end
end

println("""```
""")


# ∫ψₙ*Hψₙdx = <ψₙ|H|ψₙ> = Eₙ


println(raw"""
#### Eigen Values

```math
  \begin{aligned}
    E_n
    &=      \int \psi^\ast_n(x) \hat{H} \psi_n(x) \mathrm{d}x \\
    &=      \int \psi^\ast_n(x) \left[ \hat{V} + \hat{T} \right] \psi(x) \mathrm{d}x \\
    &=      \int \psi^\ast_n(x) \left[ V(x) - \frac{\hbar^2}{2m} \frac{\mathrm{d}^{2}}{\mathrm{d} x^{2}} \right] \psi(x) \mathrm{d}x \\
    &\simeq \int \psi^\ast_n(x) \left[ V(x)\psi(x) -\frac{\hbar^2}{2m} \frac{\psi(x+\Delta x) - 2\psi(x) + \psi(x-\Delta x)}{\Delta x^{2}} \right] \mathrm{d}x.
  \end{aligned}
```

Where, the difference formula for the 2nd-order derivative:

```math
\begin{aligned}
  % 2\psi(x)
  % + \frac{\mathrm{d}^{2} \psi(x)}{\mathrm{d} x^{2}} \Delta x^{2}
  % + O\left(\Delta x^{4}\right)
  % &=
  % \psi(x+\Delta x)
  % + \psi(x-\Delta x)
  % \\
  % \frac{\mathrm{d}^{2} \psi(x)}{\mathrm{d} x^{2}} \Delta x^{2}
  % &=
  % \psi(x+\Delta x)
  % - 2\psi(x)
  % + \psi(x-\Delta x)
  % - O\left(\Delta x^{4}\right)
  % \\
  % \frac{\mathrm{d}^{2} \psi(x)}{\mathrm{d} x^{2}}
  % &=
  % \frac{\psi(x+\Delta x) - 2\psi(x) + \psi(x-\Delta x)}{\Delta x^{2}}
  % - \frac{O\left(\Delta x^{4}\right)}{\Delta x^{2}}
  % \\
  \frac{\mathrm{d}^{2} \psi(x)}{\mathrm{d} x^{2}}
  &=
  \frac{\psi(x+\Delta x) - 2\psi(x) + \psi(x-\Delta x)}{\Delta x^{2}}
  + O\left(\Delta x^{2}\right)
\end{aligned}
```

are given by the sum of 2 Taylor series:

```math
\begin{aligned}
\psi(x+\Delta x)
&= \psi(x)
+ \frac{\mathrm{d} \psi(x)}{\mathrm{d} x} \Delta x
+ \frac{1}{2!} \frac{\mathrm{d}^{2} \psi(x)}{\mathrm{d} x^{2}} \Delta x^{2}
+ \frac{1}{3!} \frac{\mathrm{d}^{3} \psi(x)}{\mathrm{d} x^{3}} \Delta x^{3}
+ O\left(\Delta x^{4}\right),
\\
\psi(x-\Delta x)
&= \psi(x)
- \frac{\mathrm{d} \psi(x)}{\mathrm{d} x} \Delta x
+ \frac{1}{2!} \frac{\mathrm{d}^{2} \psi(x)}{\mathrm{d} x^{2}} \Delta x^{2}
- \frac{1}{3!} \frac{\mathrm{d}^{3} \psi(x)}{\mathrm{d} x^{3}} \Delta x^{3}
+ O\left(\Delta x^{4}\right).
\end{aligned}
```

```""")

@testset "∫ψₙ*Hψₙdx = <ψₙ|H|ψₙ> = Eₙ" begin
  ψHψ(HO, x; n=0, Δx=0.005) = V(HO,x)*ψ(HO,x,n=n)^2 - HO.ℏ^2/(2*HO.m)*conj(ψ(HO,x,n=n))*(ψ(HO,x+Δx,n=n)-2*ψ(HO,x,n=n)+ψ(HO,x-Δx,n=n))/Δx^2
  println("  k |  n |        analytical |         numerical ")
  println("--- | -- | ----------------- | ----------------- ")
  for k in [0.1,0.5,1.0,5.0]
  for n in 0:9
    HO = HarmonicOscillator(k=k)
    analytical = E(HO, n=n)
    numerical  = quadgk(x -> ψHψ(HO, x, n=n, Δx=0.001), -Inf, Inf, maxevals=10^3)[1]
    acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-5) : isapprox(analytical, numerical, rtol=1e-5)
    @test acceptance
    @printf("%.1f | %2d | %17.12f | %17.12f %s\n", k, n, analytical, numerical, acceptance ? "✔" : "✗")
  end
  end
end
HO = HarmonicOscillator(k=1.0, m=1.0, ℏ=1.0)

println("""```
""")
