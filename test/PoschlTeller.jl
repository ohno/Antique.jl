PT = PoschlTeller(λ=10.0)


# Pₙᵐ(x) = √(1-x²)ᵐ dᵐ/dxᵐ Pₙ(x); Pₙ(x) = 1/(2ⁿn!) dⁿ/dxⁿ (x²-1)ⁿ


println(raw"""
#### Associated Legendre Polynomials $P_n^m(x)$

```math
  \begin{aligned}
    P_n^m(x)
    &= (-1)^m \left( 1-x^2 \right)^{m/2} \frac{\mathrm{d}^m}{\mathrm{d}x^m} P_n(x) \\
    &= (-1)^m \left( 1-x^2 \right)^{m/2} \frac{\mathrm{d}^m}{\mathrm{d}x^m} \frac{1}{2^n n!} \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left[ \left( x^2-1 \right)^n \right] \\
    &= \frac{(-1)^m}{2^n} (1-x^2)^{m/2} \sum_{j=0}^{\left\lfloor\frac{n-m}{2}\right\rfloor} (-1)^j \frac{(2n-2j)!}{j! (n-j)! (n-2j-m)!} x^{(n-2j-m)}.
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
      c = (-1)^m * (1 - x^2)^(m//2) * Dm(a * Dn(b)) # Rodrigues' formula
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