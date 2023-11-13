using Antiq
using Test
using Printf
using Markdown
using QuadGK
IPW = antiq(:InfinitePotentialWell, L=1.0, m=1.0, ℏ=1.0)


# <ψᵢ|ψⱼ> = ∫ψₙ*ψₙdx = δᵢⱼ


println(raw"""
#### Normalization & Orthogonality of $\psi_n(x)$

```math
\int_{0}^{L} \psi_i^\ast(x) \psi_j(x) ~\mathrm{d}x = \delta_{ij}
```

```""")

@testset "<ψᵢ|ψⱼ> = ∫ψₙ*ψₙdx = δᵢⱼ" begin
  println("  i\t  j\tnumerical         \tanalytical        \t|error|")
  # for L in [0.1, 0.5, 1.0, 7.0]
  # for m in [0.1, 0.5, 1.0, 7.0]
  # for ℏ in [0.1, 0.5, 1.0, 7.0]
  for i in 1:10
  for j in 1:10
    numerical  = quadgk(x -> conj(IPW.ψ(x, n=i, L=IPW.L)) * IPW.ψ(x, n=j, L=IPW.L), 0.0, IPW.L, maxevals=10^3)[1]
    analytical = (i == j ? 1 : 0)
    error = analytical == 0 ? (abs(numerical) < 1e-5 ? 0.0 : Inf) : abs((numerical-analytical)/analytical)
    acceptance = error < 1e-5
    @test acceptance
    @printf("%3d\t%3d\t%.16f\t%.16f\t%.16f%%\t%s\n", i, j, numerical, analytical, error*100, acceptance ? "✔" :  "✗")
  end
  end
  # end
  # end
  # end
end

println("""```
""")


# <ψₙ|H|ψₙ>  = ∫ψₙ*Tψₙdx = Eₙ


println(raw"""
#### Eigen Values

```math
  \begin{aligned}
    E_n
    &=      \int_0^L \psi^\ast_n(x) \hat{H} \psi_n(x) ~\mathrm{d}x \\
    &=      \int_0^L \psi^\ast_n(x) \left[ \hat{V} + \hat{T} \right] \psi(x) ~\mathrm{d}x \\
    &=      \int_0^L \psi^\ast_n(x) \left[ 0 - \frac{\hbar^2}{2m} \frac{\mathrm{d}^{2}}{\mathrm{d} x^{2}} \right] \psi(x) ~\mathrm{d}x \\
    &\simeq \int_0^L \psi^\ast_n(x) \left[ -\frac{\hbar^2}{2m} \frac{\psi(x+\Delta x) - 2\psi(x) + \psi(x-\Delta x)}{\Delta x^{2}} \right] ~\mathrm{d}x.
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

@testset "<ψₙ|H|ψₙ>  = ∫ψₙ*Tψₙdx = Eₙ" begin
  ψTψ(x; n=0, L=IPW.L, m=IPW.m, ℏ=IPW.ℏ, Δx=0.01) = -ℏ^2/(2*m)*conj(IPW.ψ(x,n=n,L=L))*(IPW.ψ(x+Δx,n=n,L=L)-2*IPW.ψ(x,n=n,L=L)+IPW.ψ(x-Δx,n=n,L=L))/Δx^2
  println("  L    m    ℏ   n  numerical     \tanalytical    \t|error|")
  for L in [0.1, 1.0]
  for m in [0.1, 1.0]
  for ℏ in [0.1, 1.0]
  for n in 1:10
    numerical  = quadgk(x -> ψTψ(x, n=n, L=L, m=m, ℏ=ℏ, Δx=L*0.0001), 0, L, maxevals=10^3)[1]
    analytical = IPW.E(n=n, L=L, m=m, ℏ=ℏ)
    error = analytical == 0 ? (abs(numerical) < 1e-5 ? 0.0 : Inf) : abs((numerical-analytical)/analytical)
    acceptance = error < 1e-5
    @test acceptance
    @printf("%.1f  %.1f  %.1f  %2d  %.12f\t%.12f\t%.12f%%\t%s\n", L, m, ℏ, n, numerical, analytical, error*100, acceptance ? "✔" :  "✗")
  end
  end
  end
  end
end

println("""```""")


# <ψₙ|x|ψₙ>  = L/2


println(raw"""
#### Expected Value of $x$

```math
\langle x \rangle_{n=1}
= \int_{0}^{L} \psi_1^\ast(x) \hat{x} \psi_1(x) ~\mathrm{d}x
= \frac{2(2a)^2}{\pi^3} \left( \frac{\pi^3}{6} - \frac{\pi}{4} \right)
```
for only $n=1$.

Reference:
- [LibreTexts PHYSICS, 6.4: Expectation Values, Observables, and Uncertainty](https://phys.libretexts.org/Bookshelves/Modern_Physics/Book%3A_Spiral_Modern_Physics_(D'Alessandris)/6%3A_The_Schrodinger_Equation/6.4%3A_Expectation_Values_Observables_and_Uncertainty)
```""")

@testset "<ψₙ|x|ψₙ>  = L/2" begin
  println("   L\t  n\tnumerical         \tanalytical        \t|error|")
  for L in [0.1, 0.5, 1.0, 7.0]
  for n in 1:1
    numerical  = quadgk(x -> conj(IPW.ψ(x, n=n, L=L)) * x * IPW.ψ(x, n=n, L=L), 0, L, maxevals=10^3)[1]
    analytical = L/2
    error = analytical == 0 ? (abs(numerical) < 1e-5 ? 0.0 : Inf) : abs((numerical-analytical)/analytical)
    acceptance = error < 1e-5
    @test acceptance
    @printf("%.1f\t%3d\t%.16f\t%.16f\t%.16f%%\t%s\n", L, n, numerical, analytical, error*100, acceptance ? "✔" :  "✗")
  end
  end
end

println("""```""")


# <ψₙ|x²|ψₙ> = 2L²/π³(π³/6-π/4)


println(raw"""
#### Expected Value of $x^2$

```math
\langle x^2 \rangle_{n=1}
= \int_{0}^{L} \psi_1^\ast(x) \hat{x}^2 \psi_1(x) ~\mathrm{d}x
= \frac{2(2a)^2}{\pi^3} \left( \frac{\pi^3}{6} - \frac{\pi}{4} \right)
```

Reference:
- [LibreTexts PHYSICS, 6.4: Expectation Values, Observables, and Uncertainty](https://phys.libretexts.org/Bookshelves/Modern_Physics/Book%3A_Spiral_Modern_Physics_(D'Alessandris)/6%3A_The_Schrodinger_Equation/6.4%3A_Expectation_Values_Observables_and_Uncertainty)
```""")

@testset "<ψₙ|x²|ψₙ> = 2L²/π³(π³/6-π/4)" begin
  println("   L\t  n\tnumerical         \tanalytical        \t|error|")
  for L in [0.1, 0.5, 1.0, 7.0]
  for n in 1:1
    numerical  = quadgk(x -> conj(IPW.ψ(x, n=n, L=L)) * x^2 * IPW.ψ(x, n=n, L=L), 0, L, maxevals=10^3)[1]
    analytical = 2*L^2/π^3 * (π^3/6 - π/4)
    error = analytical == 0 ? (abs(numerical) < 1e-5 ? 0.0 : Inf) : abs((numerical-analytical)/analytical)
    acceptance = error < 1e-5
    @test acceptance
    @printf("%.1f\t%3d\t%.16f\t%.16f\t%.16f%%\t%s\n", L, n, numerical, analytical, error*100, acceptance ? "✔" :  "✗")
  end
  end
end

println("""```""")


# <ψₙ|p|ψₙ>  = ∫ψₙ*(-iℏd/dx)ψₙdx = 0


println(raw"""

#### Expected Value of $p$
```math
\langle p \rangle_{n=1}
= \int_{0}^{L} \psi_1^\ast(x) \hat{p} \psi_1(x) ~\mathrm{d}x
= 0
```

Reference:
- [LibreTexts PHYSICS, 6.4: Expectation Values, Observables, and Uncertainty](https://phys.libretexts.org/Bookshelves/Modern_Physics/Book%3A_Spiral_Modern_Physics_(D'Alessandris)/6%3A_The_Schrodinger_Equation/6.4%3A_Expectation_Values_Observables_and_Uncertainty)

---

```math
  \begin{aligned}
    \langle p \rangle_{n=1}
    &= \int_0^L \psi^\ast_n(x) \hat{p} \psi_n(x) ~\mathrm{d}x \\
    &= \int_0^L \psi^\ast_n(x) \left[ -i\hbar\frac{\mathrm{d}}{\mathrm{d} x} \right] \psi(x) ~\mathrm{d}x \\
    &\simeq \int_0^L \psi^\ast_n(x) \left[ -i\hbar \frac{\psi(x+\Delta x) - \psi(x-\Delta x)}{2\Delta x} \right] ~\mathrm{d}x.
  \end{aligned}
```

Where, the difference formula for the 2nd-order derivative:

```math
\begin{aligned}
  % 2\frac{\mathrm{d} \psi(x)}{\mathrm{d}x} \Delta x
  % + O\left(\Delta x^{3}\right)
  % &= 
  % \psi(x+\Delta x)
  % - \psi(x-\Delta x)
  % \\
  % 2\frac{\mathrm{d} \psi(x)}{\mathrm{d}x} \Delta x
  % &= 
  % \psi(x+\Delta x)
  % - \psi(x-\Delta x)
  % - O\left(\Delta x^{3}\right)
  % \\
  % \frac{\mathrm{d} \psi(x)}{\mathrm{d}x}
  % &= 
  % \frac{\psi(x+\Delta x)- \psi(x-\Delta x)}{2\Delta x}
  % - \frac{O\left(\Delta x^{3}\right)}{2\Delta x}
  % \\
  \frac{\mathrm{d} \psi(x)}{\mathrm{d}x}
  &= 
  \frac{\psi(x+\Delta x)- \psi(x-\Delta x)}{2\Delta x}
  + O\left(\Delta x^{2}\right),
\end{aligned}
```

are given by the sum of 2 Taylor series:

```math
\begin{aligned}
  \psi(x+\Delta x)
  &=
  \psi(x)
  + \frac{\mathrm{d} \psi(x)}{\mathrm{d}x} \Delta x
  + \frac{1}{2!} \frac{\mathrm{d}^{2} \psi(x)}{\mathrm{d}x^{2}} \Delta x^{2}
  + O\left(\Delta x^{3}\right),
  \\
  \psi(x-\Delta x)
  &=
  \psi(x)
  - \frac{\mathrm{d} \psi(x)}{\mathrm{d}x} \Delta x
  + \frac{1}{2!} \frac{\mathrm{d}^{2} \psi(x)}{\mathrm{d}x^{2}} \Delta x^{2}
  + O\left(\Delta x^{3}\right).
\end{aligned}
```
```""")

@testset "<ψₙ|p|ψₙ>  = ∫ψₙ*(-iℏd/dx)ψₙdx = 0" begin
  ψpψ(x; n=0, L=IPW.L, m=IPW.m, ℏ=IPW.ℏ, Δx=0.01) = -im*ℏ*conj(IPW.ψ(x,n=n,L=L))*(IPW.ψ(x+Δx,n=n,L=L)-IPW.ψ(x-Δx,n=n,L=L))/2/Δx
  println("   L\t  n\tnumerical         \tanalytical        \t|error|")
  for L in [0.1, 0.5, 1.0, 7.0]
  for n in 1:1
    numerical  = abs(quadgk(x -> ψpψ(x, n=n, L=L, m=IPW.m, ℏ=IPW.ℏ, Δx=L*0.0001), 0, L, maxevals=10^3)[1])
    analytical = 0
    error = analytical == 0 ? (abs(numerical) < 1e-5 ? 0.0 : Inf) : abs((numerical-analytical)/analytical)
    acceptance = error < 1e-5
    @test acceptance
    @printf("%.1f\t%3d\t%.16f\t%.16f\t%.16f%%\t%s\n", L, n, numerical, analytical, error*100, acceptance ? "✔" :  "✗")
  end
  end
end

println("""```""")


# <ψₙ|p²|ψₙ> = ∫ψₙ*(-ℏ²d²/dx²)ψₙdx = π²ℏ²/L²


println(raw"""
#### Expected Value of $p^2$

```math
\langle p^2 \rangle
= \int_{0}^{L} \psi_1^\ast(x) \hat{p}^2 \psi_1(x) ~\mathrm{d}x
= \frac{\pi^2\hbar^2}{L^2}
```

Reference:
- [LibreTexts PHYSICS, 6.4: Expectation Values, Observables, and Uncertainty](https://phys.libretexts.org/Bookshelves/Modern_Physics/Book%3A_Spiral_Modern_Physics_(D'Alessandris)/6%3A_The_Schrodinger_Equation/6.4%3A_Expectation_Values_Observables_and_Uncertainty)

---

```math
  \begin{aligned}
    \langle p^2 \rangle
    &= \int_0^L \psi^\ast_n(x) \hat{p} \psi_n(x) ~\mathrm{d}x \\
    &= \int_0^L \psi^\ast_n(x) \left[ -\hbar^2\frac{\mathrm{d}^2}{{\mathrm{d}x}^2} \right] \psi(x) ~\mathrm{d}x \\
    &\simeq \int_0^L \psi^\ast_n(x) \left[ -\hbar^2 \frac{\psi(x+\Delta x) - 2\psi(x) + \psi(x-\Delta x)}{\Delta x^{2}} \right] ~\mathrm{d}x.
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

@testset "<ψₙ|p²|ψₙ> = ∫ψₙ*(-ℏ²d²/dx²)ψₙdx = π²ℏ²/L²" begin
  ψp²ψ(x; n=0, L=IPW.L, m=IPW.m, ℏ=IPW.ℏ, Δx=0.01) = -ℏ^2*conj(IPW.ψ(x,n=n,L=L))*(IPW.ψ(x+Δx,n=n,L=L)-2*IPW.ψ(x,n=n,L=L)+IPW.ψ(x-Δx,n=n,L=L))/Δx^2
  println("   L\t  n\tnumerical         \tanalytical        \t|error|")
  for L in [0.1, 0.5, 1.0, 7.0]
  for n in 1:1
    numerical  = quadgk(x -> ψp²ψ(x, n=n, L=L, m=IPW.m, ℏ=IPW.ℏ, Δx=L*0.0001), 0, L, maxevals=10^3)[1]
    analytical = π^2*IPW.ℏ^2/L^2
    error = analytical == 0 ? (abs(numerical) < 1e-5 ? 0.0 : Inf) : abs((numerical-analytical)/analytical)
    acceptance = error < 1e-5
    @test acceptance
    @printf("%.1f\t%3d\t%.16f\t%.16f\t%.16f%%\t%s\n", L, n, numerical, analytical, error*100, acceptance ? "✔" :  "✗")
  end
  end
end

println("""```""")
