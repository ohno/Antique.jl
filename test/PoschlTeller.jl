PT = PoschlTeller(λ=4.0)


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
      Dn = n==0 ? x->x : Differential(x)^n          # dⁿ/dxⁿ
      Dm = m==0 ? x->x : Differential(x)^m          # dᵐ/dxᵐ
      a = 1 // (2^n * factorial(n))                 # left
      b = (x^2 - 1)^n                               # right
      c = (1 - x^2)^(m//2) * Dm(a * Dn(b)) # Rodrigues' formula
      d = expand_derivatives(c)                     # expand dⁿ/dxⁿ and dᵐ/dxᵐ
      e = simplify(d, expand=true)                  # simplify
      f = simplify(P(PT, x, n=n, m=m), expand=true) # closed-form
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


# <ψᵢ|ψⱼ> = δᵢⱼ


println(raw"""
#### Normalization & Orthogonality of $\psi_n(x)$

```math
\int \psi_i^\ast(x) \psi_j(x) \mathrm{d}x = \delta_{ij}
```

```""")

@testset "<ψᵢ|ψⱼ> = δᵢⱼ" begin
  println("  λ |   m |   ℏ |  x₀ |  i |  j |        analytical |         numerical ")
  println("--- | --- | --- | --- | -- | -- | ----------------- | ----------------- ")
  for λ in [1,2,3]
  for m in [1.0,2.0,exp(1)]
  for ℏ in [1.0,2.0,exp(1)]
  for x₀ in [1.0,2.0,exp(1)]
    PT = PoschlTeller(λ=λ, m=m, ℏ=ℏ, x₀=x₀)
    for i in 0:nₘₐₓ(PT)
    for j in 0:nₘₐₓ(PT)
      analytical = (i == j ? 1 : 0)
      numerical  = quadgk(x -> conj(ψ(PT, x, n=i)) * ψ(PT, x, n=j), -Inf, Inf, maxevals=10^3)[1]
      acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-5) : isapprox(analytical, numerical, rtol=1e-5)
      @test acceptance
      @printf("%.1f | %.1f | %.1f | %.1f | %2d | %2d | %17.12f | %17.12f %s\n", λ, m, ℏ, x₀, i, j, analytical, numerical, acceptance ? "✔" : "✗")
    end
    end
  end
  end
  end
  end
end

println("""```
""")



println(raw"""
#### Eigenvalues

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
  ψHψ(PT, x; n=0, Δx=0.005) = V(PT,x)*ψ(PT,x,n=n)^2 - PT.ℏ^2/(2*PT.m)*conj(ψ(PT,x,n=n))*(ψ(PT,x+Δx,n=n)-2*ψ(PT,x,n=n)+ψ(PT,x-Δx,n=n))/Δx^2
  println("  λ |   m |   ℏ |  x₀ |  n |        analytical |         numerical ")
  println("--- | --- | --- | --- | -- | ----------------- | ----------------- ")
  for λ in [1,2,3]
  for m in [1.0,2.0,exp(1)]
  for ℏ in [1.0,2.0,exp(1)]
  for x₀ in [1.0,2.0,exp(1)]
    PT = PoschlTeller(λ=λ, m=m, ℏ=ℏ, x₀=x₀)
    for n in 0:λ-1
      analytical = E(PT, n=n)
      numerical  = quadgk(x -> ψHψ(PT, x, n=n, Δx=0.001), -Inf, Inf, atol=1e-5)[1]
      acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-3) : isapprox(analytical, numerical, rtol=1e-5)
      @test acceptance
      @printf("%.1f | %.1f | %.1f | %.1f | %2d | %17.12f | %17.12f %s\n", λ, m, ℏ, x₀, n, analytical, numerical, acceptance ? "✔" : "✗")
    end
  end
  end
  end
  end
end
PT = PoschlTeller(λ=1.0, m=1.0, ℏ=1.0)

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
n_\mathrm{max} = \left\lfloor \lambda \right\rfloor - 1
```

```""")

@testset "0 < Eₙ₊₁ - Eₙ for 0 ≤ n ≤ nₘₐₓ" begin
  for PT in [
    PoschlTeller(λ=2, m=1.0, ℏ=1.0, x₀=1.0),
    PoschlTeller(λ=3, m=1.0, ℏ=1.0, x₀=1.0),
    PoschlTeller(λ=4, m=1.0, ℏ=1.0, x₀=1.0),
    PoschlTeller(λ=10, m=1.0, ℏ=1.0, x₀=1.0),
  ]
    @show PT
    println(" n  Eₙ          ΔE")
    for n in 0:nₘₐₓ(PT)+5
      Eₙ₊₁ = E(PT, n=n+1, nocheck=true)
      Eₙ   = E(PT, n=n, nocheck=true)
      ΔE  = Eₙ₊₁ - Eₙ
      @printf("%2d  %+9.6f  %+9.6f  ", n, Eₙ, ΔE)
      if n ≤ nₘₐₓ(PT)
        print("0 < ΔE  ")
        acceptance = 0 < ΔE
      else
        print("ΔE < 0  ")
        acceptance = ΔE < 0
      end
      @test acceptance
      println(acceptance ? "✔" : "✗")
      if n == nₘₐₓ(PT)
        println("-----------------------------------  nₘₐₓ(PT) = ", nₘₐₓ(PT))
      end
    end
    println()
  end
end

println("""```
""")
