```@meta
CurrentModule = Antique
```

# Hydrogen Atom

The hydrogen atom is the simplest 2-body Coulomb system.

## Definitions

``Z`` is the atomic number. The domains of the potential and the wave functions are $0\leq r \lt \infty, 0\leq \theta \lt \pi, 0\leq \varphi \lt 2\pi$.

#### Schrödinger Equation
```math
  \hat{H} \psi(\pmb{r}) = E \psi(\pmb{r})
```

#### Hamiltonian
```math
  \hat{H} = - \frac{\hbar^2}{2\mu} \frac{\mathrm{d}^2}{\mathrm{d}r ^2} + V(r).
```
Where, $\mu=\left(\frac{1}{m_\mathrm{e}}+\frac{1}{m_\mathrm{p}}\right)^{-1}$ is the reduced mass of electron $\mathrm{e}$ and proton $\mathrm{p}$. $\mu = m_\mathrm{e}$ holds in the limit $m_\mathrm{p}\rightarrow\infty$. 

#### Potential
`V(r; Z=Z, a₀=a₀)`
```math
  \begin{aligned}
  V(r)
  &= - \frac{Ze^2}{4\pi\varepsilon_0 r} 
  &= - \frac{e^2}{4\pi\varepsilon_0 a_0} \frac{Z}{r/a_0}
  &= - \frac{Z}{r/a_0} E_\mathrm{h}
  \end{aligned}
```

#### Eigen Values
`E(; n=1, Z=Z, Eₕ=Eₕ)`
```math
  E_n = -\frac{m_\mathrm{e} e^4 Z^2}{2n^2(4\pi\varepsilon_0)^2\hbar^2} = -\frac{Z^2}{2n^2} E_\mathrm{h}
```
Where, ``E_\mathrm{h}`` is the Hartree energy, one of atomic unit. About atomic units, see section 3.9.2 of the [IUPAC GreenBook](https://iupac.org/what-we-do/books/greenbook/). In other units, ``E_\mathrm{h} = 27.211~386~245~988(53)~\mathrm{eV}`` from [here](https://physics.nist.gov/cgi-bin/cuu/Value?hrev).

#### Eigen Functions
`ψ(r, θ, φ; n=1, l=0, m=0, Z=Z, a₀=a₀)`
```math
  \psi_{nlm}(\pmb{r}) = R_{nl}(r) Y_{lm}(\theta,\varphi)
```

#### Radial Functions
`R(r; n=1, l=0, Z=Z, a₀=a₀)`
```math
    R_{nl}(r) = -\sqrt{\frac{(n-l-1)!}{2n(n+l)!} \left(\frac{2Z}{n a_0}\right)^3} \left(\frac{2Zr}{n a_0}\right)^l \exp \left(-\frac{Zr}{n a_0}\right) L_{n+l}^{2l+1} \left(\frac{2Zr}{n a_0}\right)
```
Where, Laguerre polynomials are defined as ``L_n(x) = \frac{1}{n!} \mathrm{e}^x \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left( \mathrm{e}^{-x} x^n \right)``, and associated Laguerre polynomials are defined as ``L_n^{k}(x) = \frac{\mathrm{d}^k}{\mathrm{d}x^k} L_n(x)``. Note that, replace ``2n(n+l)!`` with ``2n[(n+l)!]^3`` if Laguerre polynomials are defined as ``L_n(x) = \mathrm{e}^x \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left( \mathrm{e}^{-x} x^n \right)``.

#### Associated Laguerre Polynomials
`L(x; n=0, k=0)`

Rodrigues' formula & closed-form:
```math
  \begin{aligned}
  L_n^{k}(x)
    &= \frac{\mathrm{d}^k}{\mathrm{d}x^k} L_n(x) \\
    &= \frac{\mathrm{d}^k}{\mathrm{d}x^k} \frac{1}{n!} \mathrm{e}^x \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left( \mathrm{e}^{-x} x^n \right) \\
    &= \sum_{m=0}^{n-k} (-1)^{m+k} \frac{n!}{m!(m+k)!(n-m-k)!} x^m \\
    &= (-1)^k L_{n-k}^{(k)}(x)
  \end{aligned}
```

Where, Laguerre polynomials are defined as ``L_n(x)=\frac{1}{n!}\mathrm{e}^x \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left( \mathrm{e}^{-x} x^n \right)``.

Examples:
```math
  \begin{aligned}
    L_0^0(x) &= 1, \\
    L_1^0(x) &= 1 - x, \\
    L_1^1(x) &= 1, \\
    L_2^0(x) &= 1 - 2 x + 1/2 x^2, \\
    L_2^1(x) &= 2 - x, \\
    L_2^2(x) &= 1, \\
    L_3^0(x) &= 1 - 3 x + 3/2 x^2 - 1/6 x^3, \\
    L_3^1(x) &= 3 - 3 x + 1/2 x^2, \\
    L_3^2(x) &= 3 - x, \\
    L_3^3(x) &= 1, \\
    L_4^0(x) &= 1 - 4 x + 3 x^2 - 2/3 x^3 + 5/12 x^4, \\
    L_4^1(x) &= 4 - 6 x + 2 x^2 - 1/6 x^3, \\
    L_4^2(x) &= 6 - 4 x + 1/2 x^2, \\
    L_4^3(x) &= 4 - x, \\
    L_4^4(x) &= 1, \\
    \vdots.
  \end{aligned}
```

#### Spherical Harmonics
`Y(θ, φ; l=0, m=0)`
```math
    Y_{lm}(\theta,\varphi) = (-1)^{\frac{|m|+m}{2}} \sqrt{\frac{2l+1}{4\pi} \frac{(l-|m|)!}{(l+|m|)!}} P_l^{|m|} (\cos\theta) \mathrm{e}^{im\varphi}
```
Note that some variants are connected by 
```math
i^{|m|+m} \sqrt{\frac{(l-|m|)!}{(l+|m|)!}} P_l^{|m|} = (-1)^{\frac{|m|+m}{2}} \sqrt{\frac{(l-|m|)!}{(l+|m|)!}} P_l^{|m|} = (-1)^m \sqrt{\frac{(l-m)!}{(l+m)!}} P_l^{m}.
```

#### Associated Legendre Polynomials
`P(x; n=0, m=0)`

Rodrigues' formula & closed-form:
```math
  \begin{aligned}
    P_n^m(x)
    &= \left( 1-x^2 \right)^{m/2} \frac{\mathrm{d}^m}{\mathrm{d}x^m} P_n(x) \\
    &= \left( 1-x^2 \right)^{m/2} \frac{\mathrm{d}^m}{\mathrm{d}x^m} \frac{1}{2^n n!} \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left[ \left( x^2-1 \right)^n \right] \\
    &= \frac{1}{2^n} (1-x^2)^{m/2} \sum_{j=0}^{\left\lfloor\frac{n-m}{2}\right\rfloor} (-1)^j \frac{(2n-2j)!}{j! (n-j)! (n-2j-m)!} x^{(n-2j-m)}.
  \end{aligned}
```

Where, Legendre polynomials are defined as ``P_n(x) = \frac{1}{2^n n!} \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left[ \left( x^2-1 \right)^n \right]``. Note that ``P_l^{-m} = (-1)^m \frac{(l-m)!}{(l+m)!} P_l^m`` for ``m<0``. (It is not compatible with ``P_k^m(t) = (-1)^m\left( 1-t^2 \right)^{m/2} \frac{\mathrm{d}^m P_k(t)}{\mathrm{d}t^m}`` caused by ``(-1)^m``.) The specific formulae are given below.

Examples:
```math
  \begin{aligned}
    P_{0}^{0}(x) &= 1, \\
    P_{1}^{0}(x) &= x, \\
    P_{1}^{1}(x) &= \left(+1\right)\sqrt{1-x^2}, \\
    P_{2}^{0}(x) &= -1/2 + 3/2 x^{2}, \\
    P_{2}^{1}(x) &= \left(-3 x\right)\sqrt{1-x^2}, \\
    P_{2}^{2}(x) &= 3 - 6 x, \\
    P_{3}^{0}(x) &= -3/2 x + 5/2 x^{3}, \\
    P_{3}^{1}(x) &= \left(3/2 - 15/2 x^{2}\right)\sqrt{1-x^2}, \\
    P_{3}^{2}(x) &= 15 x - 30 x^{2}, \\
    P_{3}^{3}(x) &= \left(15 - 30 x\right)\sqrt{1-x^2}, \\
    P_{4}^{0}(x) &= 3/8 - 15/4 x^{2} + 35/8 x^{4}, \\
    P_{4}^{1}(x) &= \left(- 15/2 x + 35/2 x^{3}\right)\sqrt{1-x^2}, \\
    P_{4}^{2}(x) &= -15/2 + 15 x + 105/2 x^{2} - 105 x^{3}, \\
    P_{4}^{3}(x) &= \left(105 x - 210 x^{2}\right)\sqrt{1-x^2}, \\
    P_{4}^{4}(x) &= 105 - 420 x + 420 x^{2}, \\
    & \vdots.
  \end{aligned}
```

#### References
- cpprefjp, [legendre](https://cpprefjp.github.io/reference/cmath/legendre.html), [assoc_legendre](https://cpprefjp.github.io/reference/cmath/assoc_legendre.html), [laguerre](https://cpprefjp.github.io/reference/cmath/laguerre.html), [assoc_laguerre](https://cpprefjp.github.io/reference/cmath/assoc_laguerre.html)
- The Digital Library of Mathematical Functions (DLMF), [18.3 Table1](https://dlmf.nist.gov/18.3#T1), [18.5 Table1](https://dlmf.nist.gov/18.5#T1), [18.5.16](https://dlmf.nist.gov/18.5#E16), [18.3 Table1](https://dlmf.nist.gov/18.3#T1), [18.5 Table1](https://dlmf.nist.gov/18.5#T1), [18.5.17](https://dlmf.nist.gov/18.5#E17), [18.3 Table1](https://dlmf.nist.gov/18.3#T1), [18.5 Table1](https://dlmf.nist.gov/18.5#T1), [18.5.12](https://dlmf.nist.gov/18.5#E12)
- L. D. Landau, E. M. Lifshitz, Quantum Mechanics (Pergamon Press, 1965), [p.598 (c.1)](https://archive.org/details/ost-physics-landaulifshitz-quantummechanics/page/n611/mode/2up), [p.598 (c.4)](https://archive.org/details/ost-physics-landaulifshitz-quantummechanics/page/n611/mode/2up), [p.603 (d.13)](https://archive.org/details/ost-physics-landaulifshitz-quantummechanics/page/n615/mode/2up), [p.603 (d.13)](https://archive.org/details/ost-physics-landaulifshitz-quantummechanics/page/n615/mode/2up)
- L. I. Schiff, Quantum Mechanics (McGraw-Hill Book Company, 1968), [p.79 (14.12)](https://archive.org/details/ost-physics-schiff-quantummechanics/page/n95/mode/1up), [p.93 (16.19)](https://archive.org/details/ost-physics-schiff-quantummechanics/page/n109/mode/1up)
- A. Messiah, Quanfum Mechanics (Dover Publications, 1999), [p.493 (B.72)](https://archive.org/details/quantummechanics0000mess/page/491/mode/1up), [p.494 Table](https://archive.org/details/quantummechanics0000mess/page/494/mode/1up), [p.493 (B.72)](https://archive.org/details/quantummechanics0000mess/page/491/mode/1up), [p.483 (B.12)](https://archive.org/details/quantummechanics0000mess/page/483/mode/1up), [p.483 (B.12)](https://archive.org/details/quantummechanics0000mess/page/483/mode/1up)
- W. Greiner, Quantum Mechanics: An Introduction Third Edition (Springer, 1994), [p.83 (4)](https://archive.org/details/quantummechanics0001grei_u4x0/page/83/mode/1up), [p.83 (5)](https://archive.org/details/quantummechanics0001grei_u4x0/page/83/mode/1up), [p.149 (21)](https://archive.org/details/quantummechanics0001grei_u4x0/page/149/mode/1up)
- D. J. Griffiths, Introduction to Quantum Mechanics (Prentice Hall, 1995), [p.126 (4.28)](https://archive.org/details/griffiths-introduction-to-quantum-mechanics/page/126/mode/1up), [p.96 Table3.1](https://archive.org/details/griffiths-introduction-to-quantum-mechanics/page/95/mode/1up), [p.126 (4.27)](https://archive.org/details/griffiths-introduction-to-quantum-mechanics/page/126/mode/1up), [p.139 (4.88)](https://archive.org/details/griffiths-introduction-to-quantum-mechanics/page/139/mode/1up), [p.140 Table4.4](https://archive.org/details/griffiths-introduction-to-quantum-mechanics/page/140/mode/1up), [p.139 (4.87)](https://archive.org/details/griffiths-introduction-to-quantum-mechanics/page/139/mode/1up), [p.140 Table4.5](https://archive.org/details/griffiths-introduction-to-quantum-mechanics/page/140/mode/1up)
- D. A. McQuarrie, J. D. Simon, Physical Chemistry: A Molecular Approach (University Science Books, 1997), [p.195 Table6.1](https://archive.org/details/McQuarrieSimonPhysicalChemistrySolutions/McQuarrie_Simon_Physical_Chemistry1997/page/n218/mode/1up), [p.196 (6.26)](https://archive.org/details/McQuarrieSimonPhysicalChemistrySolutions/McQuarrie_Simon_Physical_Chemistry1997/page/n219/mode/1up), [p.196 Table6.2](https://archive.org/details/McQuarrieSimonPhysicalChemistrySolutions/McQuarrie_Simon_Physical_Chemistry1997/page/n220/mode/1up), [p.207 Table6.4](https://archive.org/details/McQuarrieSimonPhysicalChemistrySolutions/McQuarrie_Simon_Physical_Chemistry1997/page/n230/mode/1up)
- P. W. Atkins, J. De Paula, Atkins' Physical Chemistry, 8th edition (W. H. Freeman, 2008), [p.234](https://archive.org/details/atkinsphysicalch00pwat/page/324/mode/2up?q=Laguerre)
- [J. J. Sakurai, J. Napolitano, Modern Quantum Mechanics Third Edition (Cambridge University Press, 2021)](https://doi.org/10.1017/9781108587280), p.245 Problem 3.30.b, 

## Usage & Examples

[Install Antique.jl](@ref Install) for the first run and run `using Antique` before each use. The function `antique(model, parameters...)` returns a module that has `E()`, `ψ(r)`, `V(r)` and some other functions. In this system, the model name is specified by `:HydrogenAtom` and several parameters `Z`, `Eₕ`, `mₑ`, `a₀` and `ℏ` are set as optional arguments.

```julia
using Antique
H = antique(:HydrogenAtom, Z=1, Eₕ=1.0, a₀=1.0, mₑ=1.0, ℏ=1.0)
```




Parameters:

```julia
julia> H.Z
1

julia> H.Eₕ
1.0

julia> H.mₑ
1.0

julia> H.a₀
1.0

julia> H.ℏ
1.0
```



Eigen values:

```julia
julia> H.E(n=1)
-0.5

julia> H.E(n=2)
-0.125
```



Wave length ($n=2\rightarrow1$, the first line of the Lyman series):

```julia
Eₕ2nm⁻¹ = 2.1947463136320e-2 # https://physics.nist.gov/cgi-bin/cuu/CCValue?hrminv
println("ΔE = ", H.E(n=2) - H.E(n=1), " Eₕ")
println("λ  = ", ((H.E(n=2)-H.E(n=1))*Eₕ2nm⁻¹)^-1, " nm")
```

```
ΔE = 0.375 Eₕ
λ  = 121.50227341098497 nm
```





Hyperfine Splitting:

```julia
# constants: https://doi.org/10.1103/RevModPhys.93.025010
e  = 1.602176634e-19    # C      https://physics.nist.gov/cgi-bin/cuu/Value?e
h  = 6.62607015e-34     # J Hz-1 https://physics.nist.gov/cgi-bin/cuu/Value?h
c  = 299792458          # m s-1  https://physics.nist.gov/cgi-bin/cuu/Value?c
a0 = 5.29177210903e-11  # m      https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0
µ0 = 1.25663706212e-6   # N A-2  https://physics.nist.gov/cgi-bin/cuu/Value?mu0
µB = 9.2740100783e-24   # J T-1  https://physics.nist.gov/cgi-bin/cuu/Value?mub
µN = 5.0507837461e-27   # J T-1  https://physics.nist.gov/cgi-bin/cuu/Value?mun
ge = 2.00231930436256   #        https://physics.nist.gov/cgi-bin/cuu/Value?gem
gp = 5.5856946893       #        https://physics.nist.gov/cgi-bin/cuu/Value?gp

# calculation: https://doi.org/10.1119/1.12733
δ = abs(H.ψ(0,0,0))^2
ΔE = 2 / 3 * µ0 * µN * µB * gp * ge * δ * a0^(-3)
println("1/π    = ", 1/π)
println("<δ(r)> = ", δ, " a₀⁻³")
println("<δ(r)> = ", δ * a0^(-3), " m⁻³")
println("ΔE = ", ΔE, " J")
println("ν = ΔE/h = ", ΔE / h * 1e-6, " MHz")
println("λ = hc/ΔE = ", h*c/ΔE*100, " cm")
```

```
1/π    = 0.3183098861837907
<δ(r)> = 0.3183098861837908 a₀⁻³
<δ(r)> = 2.1480615849063944e30 m⁻³
ΔE = 9.427622831641132e-25 J
ν = ΔE/h = 1422.8075794882932 MHz
λ = hc/ΔE = 21.070485027063118 cm
```





Potential energy curve:

```julia
using Plots
plot(xlims=(0.0,15.0), ylims=(-0.6,0.05), xlabel="\$r~/~a_0\$", ylabel="\$V(r)/E_\\mathrm{h},~E_n/E_\\mathrm{h}\$", legend=:bottomright, size=(480,400))
plot!(0.1:0.01:15, r -> H.V(r), lc=:black, lw=2, label="") # potential
```

![](./assets/fig//HydrogenAtom_6_1.png)



Potential energy curve, Energy levels:

```julia
using Plots
plot(xlims=(0.0,15.0), ylims=(-0.6,0.05), xlabel="\$r~/~a_0\$", ylabel="\$V(r)/E_\\mathrm{h}\$", legend=:bottomright, size=(480,400))
for n in 0:10
  plot!(0.0:0.01:15, r -> H.E(n=n) > H.V(r) ? H.E(n=n) : NaN, lc=n, lw=1, label="") # energy level
end
plot!(0.1:0.01:15, r -> H.V(r), lc=:black, lw=2, label="") # potential
```

![](./assets/fig//HydrogenAtom_7_1.png)



Radial functions:

```julia
using Plots
plot(xlabel="\$r~/~a_0\$", ylabel="\$r^2|R_{nl}(r)|^2~/~a_0^{-1}\$", ylims=(-0.01,0.55), xticks=0:1:20, size=(480,400), dpi=300)
for n in 1:3
  for l in 0:n-1
    plot!(0:0.01:20, r->r^2*H.R(r,n=n,l=l)^2, lc=n, lw=2, ls=[:solid,:dash,:dot,:dashdot,:dashdotdot][l+1], label="\$n = $n, l=$l\$")
  end
end
plot!()
```

![](./assets/fig//HydrogenAtom_8_1.png)



## Testing

Unit testing and Integration testing were done using computer algebra system ([Symbolics.jl](https://symbolics.juliasymbolics.org/stable/)) and numerical integration ([QuadGK.jl](https://juliamath.github.io/QuadGK.jl/stable/)). The test script is [here](https://github.com/ohno/Antique.jl/blob/main/test/HydrogenAtom.jl).

#### Associated Legendre Polynomials $P_n^m(x)$

```math
  \begin{aligned}
    P_n^m(x)
    &= \left( 1-x^2 \right)^{m/2} \frac{\mathrm{d}^m}{\mathrm{d}x^m} P_n(x) \\
    &= \left( 1-x^2 \right)^{m/2} \frac{\mathrm{d}^m}{\mathrm{d}x^m} \frac{1}{2^n n!} \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left[ \left( x^2-1 \right)^n \right] \\
    &= \frac{1}{2^n} (1-x^2)^{m/2} \sum_{j=0}^{\left\lfloor\frac{n-m}{2}\right\rfloor} (-1)^j \frac{(2n-2j)!}{j! (n-j)! (n-2j-m)!} x^{(n-2j-m)}.
  \end{aligned}
```

``n=0, m=0:`` ✔
```math
\begin{aligned}
  P_{0}^{0}(x)
    = 1
  &= 1 \\
  &= 1
\end{aligned}
```

``n=1, m=0:`` ✔
```math
\begin{aligned}
  P_{1}^{0}(x)
    = \frac{1}{2} \frac{\mathrm{d}}{\mathrm{d}x} \left( -1 + x^{2} \right)
  &= x \\
  &= x
\end{aligned}
```

``n=1, m=1:`` ✔
```math
\begin{aligned}
  P_{1}^{1}(x)
    = \left( 1 - x^{2} \right)^{\frac{1}{2}} \frac{\mathrm{d}}{\mathrm{d}x} \frac{1}{2} \frac{\mathrm{d}}{\mathrm{d}x} \left( -1 + x^{2} \right)
  &= \left( 1 - x^{2} \right)^{\frac{1}{2}} \\
  &= \left( 1 - x^{2} \right)^{\frac{1}{2}}
\end{aligned}
```

``n=2, m=0:`` ✔
```math
\begin{aligned}
  P_{2}^{0}(x)
    = \frac{1}{8} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \left( -1 + x^{2} \right)^{2}
  &= \frac{-1}{2} + \frac{3}{2} x^{2} \\
  &= \frac{-1}{2} + \frac{3}{2} x^{2}
\end{aligned}
```

``n=2, m=1:`` ✔
```math
\begin{aligned}
  P_{2}^{1}(x)
    = \left( 1 - x^{2} \right)^{\frac{1}{2}} \frac{\mathrm{d}}{\mathrm{d}x} \frac{1}{8} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \left( -1 + x^{2} \right)^{2}
  &= 3 \left( 1 - x^{2} \right)^{\frac{1}{2}} x \\
  &= 3 \left( 1 - x^{2} \right)^{\frac{1}{2}} x
\end{aligned}
```

``n=2, m=2:`` ✔
```math
\begin{aligned}
  P_{2}^{2}(x)
    = \left( 1 - x^{2} \right) \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{1}{8} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \left( -1 + x^{2} \right)^{2}
  &= 3 - 3 x^{2} \\
  &= 3 - 3 x^{2}
\end{aligned}
```

``n=3, m=0:`` ✔
```math
\begin{aligned}
  P_{3}^{0}(x)
    = \frac{1}{48} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \left( -1 + x^{2} \right)^{3}
  &=  - \frac{3}{2} x + \frac{5}{2} x^{3} \\
  &=  - \frac{3}{2} x + \frac{5}{2} x^{3}
\end{aligned}
```

``n=3, m=1:`` ✔
```math
\begin{aligned}
  P_{3}^{1}(x)
    = \left( 1 - x^{2} \right)^{\frac{1}{2}} \frac{\mathrm{d}}{\mathrm{d}x} \frac{1}{48} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \left( -1 + x^{2} \right)^{3}
  &=  - \frac{3}{2} \left( 1 - x^{2} \right)^{\frac{1}{2}} + \frac{15}{2} \left( 1 - x^{2} \right)^{\frac{1}{2}} x^{2} \\
  &=  - \frac{3}{2} \left( 1 - x^{2} \right)^{\frac{1}{2}} + \frac{15}{2} \left( 1 - x^{2} \right)^{\frac{1}{2}} x^{2}
\end{aligned}
```

``n=3, m=2:`` ✔
```math
\begin{aligned}
  P_{3}^{2}(x)
    = \left( 1 - x^{2} \right) \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{1}{48} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \left( -1 + x^{2} \right)^{3}
  &= 15 x - 15 x^{3} \\
  &= 15 x - 15 x^{3}
\end{aligned}
```

``n=3, m=3:`` ✔
```math
\begin{aligned}
  P_{3}^{3}(x)
    = \left( 1 - x^{2} \right)^{\frac{3}{2}} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{1}{48} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \left( -1 + x^{2} \right)^{3}
  &= 15 \left( 1 - x^{2} \right)^{\frac{3}{2}} \\
  &= 15 \left( 1 - x^{2} \right)^{\frac{3}{2}}
\end{aligned}
```

``n=4, m=0:`` ✔
```math
\begin{aligned}
  P_{4}^{0}(x)
    = \frac{1}{384} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \left( -1 + x^{2} \right)^{4}
  &= \frac{3}{8} - \frac{15}{4} x^{2} + \frac{35}{8} x^{4} \\
  &= \frac{3}{8} - \frac{15}{4} x^{2} + \frac{35}{8} x^{4}
\end{aligned}
```

``n=4, m=1:`` ✔
```math
\begin{aligned}
  P_{4}^{1}(x)
    = \left( 1 - x^{2} \right)^{\frac{1}{2}} \frac{\mathrm{d}}{\mathrm{d}x} \frac{1}{384} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \left( -1 + x^{2} \right)^{4}
  &=  - \frac{15}{2} \left( 1 - x^{2} \right)^{\frac{1}{2}} x + \frac{35}{2} \left( 1 - x^{2} \right)^{\frac{1}{2}} x^{3} \\
  &=  - \frac{15}{2} \left( 1 - x^{2} \right)^{\frac{1}{2}} x + \frac{35}{2} \left( 1 - x^{2} \right)^{\frac{1}{2}} x^{3}
\end{aligned}
```

``n=4, m=2:`` ✔
```math
\begin{aligned}
  P_{4}^{2}(x)
    = \left( 1 - x^{2} \right) \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{1}{384} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \left( -1 + x^{2} \right)^{4}
  &= \frac{-15}{2} + 60 x^{2} - \frac{105}{2} x^{4} \\
  &= \frac{-15}{2} + 60 x^{2} - \frac{105}{2} x^{4}
\end{aligned}
```

``n=4, m=3:`` ✔
```math
\begin{aligned}
  P_{4}^{3}(x)
    = \left( 1 - x^{2} \right)^{\frac{3}{2}} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{1}{384} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \left( -1 + x^{2} \right)^{4}
  &= 105 \left( 1 - x^{2} \right)^{\frac{3}{2}} x \\
  &= 105 \left( 1 - x^{2} \right)^{\frac{3}{2}} x
\end{aligned}
```

``n=4, m=4:`` ✔
```math
\begin{aligned}
  P_{4}^{4}(x)
    = \left( 1 - x^{2} \right)^{2} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{1}{384} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \left( -1 + x^{2} \right)^{4}
  &= 105 \left( 1 - x^{2} \right)^{2} \\
  &= 105 \left( 1 - x^{2} \right)^{2}
\end{aligned}
```

```
Test Summary:                                                   | Pass  Total  Time
Pₙᵐ(x) = √(1-x²)ᵐ dᵐ/dxᵐ Pₙ(x); Pₙ(x) = 1/(2ⁿn!) dⁿ/dxⁿ (x²-1)ⁿ |   15     15  1.2s
```

#### Normalization & Orthogonality of $P_n^m(x)$

```math
\int_{-1}^{1} P_i^m(x) P_j^m(x) \mathrm{d}x = \frac{2(j+m)!}{(2j+1)(j-m)!} \delta_{ij}
```
```
 m |  i |  j |        analytical |         numerical 
-- | -- | -- | ----------------- | ----------------- 
 0 |  0 |  0 |    2.000000000000 |    2.000000000000 ✔
 0 |  0 |  1 |    0.000000000000 |    0.000000000000 ✔
 0 |  0 |  2 |    0.000000000000 |    0.000000000000 ✔
 0 |  0 |  3 |    0.000000000000 |   -0.000000000000 ✔
 0 |  0 |  4 |    0.000000000000 |    0.000000000000 ✔
 0 |  0 |  5 |    0.000000000000 |    0.000000000000 ✔
 0 |  0 |  6 |    0.000000000000 |    0.000000000000 ✔
 0 |  0 |  7 |    0.000000000000 |    0.000000000000 ✔
 0 |  0 |  8 |    0.000000000000 |    0.000000000000 ✔
 0 |  0 |  9 |    0.000000000000 |   -0.000000000000 ✔
 0 |  1 |  0 |    0.000000000000 |    0.000000000000 ✔
 0 |  1 |  1 |    0.666666666667 |    0.666666666667 ✔
 0 |  1 |  2 |    0.000000000000 |    0.000000000000 ✔
 0 |  1 |  3 |    0.000000000000 |   -0.000000000000 ✔
 0 |  1 |  4 |    0.000000000000 |    0.000000000000 ✔
 0 |  1 |  5 |    0.000000000000 |    0.000000000000 ✔
 0 |  1 |  6 |    0.000000000000 |    0.000000000000 ✔
 0 |  1 |  7 |    0.000000000000 |    0.000000000000 ✔
 0 |  1 |  8 |    0.000000000000 |    0.000000000000 ✔
 0 |  1 |  9 |    0.000000000000 |   -0.000000000000 ✔
 0 |  2 |  0 |    0.000000000000 |    0.000000000000 ✔
 0 |  2 |  1 |    0.000000000000 |    0.000000000000 ✔
 0 |  2 |  2 |    0.400000000000 |    0.400000000000 ✔
 0 |  2 |  3 |    0.000000000000 |    0.000000000000 ✔
 0 |  2 |  4 |    0.000000000000 |    0.000000000000 ✔
 0 |  2 |  5 |    0.000000000000 |    0.000000000000 ✔
 0 |  2 |  6 |    0.000000000000 |    0.000000000000 ✔
 0 |  2 |  7 |    0.000000000000 |    0.000000000000 ✔
 0 |  2 |  8 |    0.000000000000 |    0.000000000000 ✔
 0 |  2 |  9 |    0.000000000000 |   -0.000000000000 ✔
 0 |  3 |  0 |    0.000000000000 |   -0.000000000000 ✔
 0 |  3 |  1 |    0.000000000000 |   -0.000000000000 ✔
 0 |  3 |  2 |    0.000000000000 |    0.000000000000 ✔
 0 |  3 |  3 |    0.285714285714 |    0.285714285714 ✔
 0 |  3 |  4 |    0.000000000000 |    0.000000000000 ✔
 0 |  3 |  5 |    0.000000000000 |   -0.000000000000 ✔
 0 |  3 |  6 |    0.000000000000 |    0.000000000000 ✔
 0 |  3 |  7 |    0.000000000000 |    0.000000000000 ✔
 0 |  3 |  8 |    0.000000000000 |   -0.000000000000 ✔
 0 |  3 |  9 |    0.000000000000 |   -0.000000000000 ✔
 0 |  4 |  0 |    0.000000000000 |    0.000000000000 ✔
 0 |  4 |  1 |    0.000000000000 |    0.000000000000 ✔
 0 |  4 |  2 |    0.000000000000 |    0.000000000000 ✔
 0 |  4 |  3 |    0.000000000000 |    0.000000000000 ✔
 0 |  4 |  4 |    0.222222222222 |    0.222222222222 ✔
 0 |  4 |  5 |    0.000000000000 |   -0.000000000000 ✔
 0 |  4 |  6 |    0.000000000000 |   -0.000000000000 ✔
 0 |  4 |  7 |    0.000000000000 |   -0.000000000000 ✔
 0 |  4 |  8 |    0.000000000000 |    0.000000000000 ✔
 0 |  4 |  9 |    0.000000000000 |    0.000000000000 ✔
 0 |  5 |  0 |    0.000000000000 |    0.000000000000 ✔
 0 |  5 |  1 |    0.000000000000 |    0.000000000000 ✔
 0 |  5 |  2 |    0.000000000000 |    0.000000000000 ✔
 0 |  5 |  3 |    0.000000000000 |   -0.000000000000 ✔
 0 |  5 |  4 |    0.000000000000 |   -0.000000000000 ✔
 0 |  5 |  5 |    0.181818181818 |    0.181818181818 ✔
 0 |  5 |  6 |    0.000000000000 |    0.000000000000 ✔
 0 |  5 |  7 |    0.000000000000 |    0.000000000000 ✔
 0 |  5 |  8 |    0.000000000000 |   -0.000000000000 ✔
 0 |  5 |  9 |    0.000000000000 |   -0.000000000000 ✔
 0 |  6 |  0 |    0.000000000000 |    0.000000000000 ✔
 0 |  6 |  1 |    0.000000000000 |    0.000000000000 ✔
 0 |  6 |  2 |    0.000000000000 |    0.000000000000 ✔
 0 |  6 |  3 |    0.000000000000 |    0.000000000000 ✔
 0 |  6 |  4 |    0.000000000000 |   -0.000000000000 ✔
 0 |  6 |  5 |    0.000000000000 |    0.000000000000 ✔
 0 |  6 |  6 |    0.153846153846 |    0.153846153846 ✔
 0 |  6 |  7 |    0.000000000000 |   -0.000000000000 ✔
 0 |  6 |  8 |    0.000000000000 |    0.000000000000 ✔
 0 |  6 |  9 |    0.000000000000 |    0.000000000000 ✔
 0 |  7 |  0 |    0.000000000000 |    0.000000000000 ✔
 0 |  7 |  1 |    0.000000000000 |    0.000000000000 ✔
 0 |  7 |  2 |    0.000000000000 |    0.000000000000 ✔
 0 |  7 |  3 |    0.000000000000 |    0.000000000000 ✔
 0 |  7 |  4 |    0.000000000000 |   -0.000000000000 ✔
 0 |  7 |  5 |    0.000000000000 |    0.000000000000 ✔
 0 |  7 |  6 |    0.000000000000 |   -0.000000000000 ✔
 0 |  7 |  7 |    0.133333333333 |    0.133333333333 ✔
 0 |  7 |  8 |    0.000000000000 |    0.000000000000 ✔
 0 |  7 |  9 |    0.000000000000 |   -0.000000000000 ✔
 0 |  8 |  0 |    0.000000000000 |    0.000000000000 ✔
 0 |  8 |  1 |    0.000000000000 |    0.000000000000 ✔
 0 |  8 |  2 |    0.000000000000 |    0.000000000000 ✔
 0 |  8 |  3 |    0.000000000000 |   -0.000000000000 ✔
 0 |  8 |  4 |    0.000000000000 |    0.000000000000 ✔
 0 |  8 |  5 |    0.000000000000 |   -0.000000000000 ✔
 0 |  8 |  6 |    0.000000000000 |    0.000000000000 ✔
 0 |  8 |  7 |    0.000000000000 |    0.000000000000 ✔
 0 |  8 |  8 |    0.117647058824 |    0.117647058824 ✔
 0 |  8 |  9 |    0.000000000000 |    0.000000000000 ✔
 0 |  9 |  0 |    0.000000000000 |   -0.000000000000 ✔
 0 |  9 |  1 |    0.000000000000 |   -0.000000000000 ✔
 0 |  9 |  2 |    0.000000000000 |   -0.000000000000 ✔
 0 |  9 |  3 |    0.000000000000 |   -0.000000000000 ✔
 0 |  9 |  4 |    0.000000000000 |    0.000000000000 ✔
 0 |  9 |  5 |    0.000000000000 |   -0.000000000000 ✔
 0 |  9 |  6 |    0.000000000000 |    0.000000000000 ✔
 0 |  9 |  7 |    0.000000000000 |   -0.000000000000 ✔
 0 |  9 |  8 |    0.000000000000 |    0.000000000000 ✔
 0 |  9 |  9 |    0.105263157895 |    0.105263157895 ✔
 1 |  1 |  1 |    1.333333333333 |    1.333333333333 ✔
 1 |  1 |  2 |    0.000000000000 |    0.000000000000 ✔
 1 |  1 |  3 |    0.000000000000 |    0.000000000000 ✔
 1 |  1 |  4 |    0.000000000000 |    0.000000000000 ✔
 1 |  1 |  5 |    0.000000000000 |    0.000000000000 ✔
 1 |  1 |  6 |    0.000000000000 |   -0.000000000000 ✔
 1 |  1 |  7 |    0.000000000000 |    0.000000000000 ✔
 1 |  1 |  8 |    0.000000000000 |    0.000000000000 ✔
 1 |  1 |  9 |    0.000000000000 |    0.000000000000 ✔
 1 |  2 |  1 |    0.000000000000 |    0.000000000000 ✔
 1 |  2 |  2 |    2.400000000000 |    2.400000000000 ✔
 1 |  2 |  3 |    0.000000000000 |    0.000000000000 ✔
 1 |  2 |  4 |    0.000000000000 |    0.000000000000 ✔
 1 |  2 |  5 |    0.000000000000 |   -0.000000000000 ✔
 1 |  2 |  6 |    0.000000000000 |    0.000000000000 ✔
 1 |  2 |  7 |    0.000000000000 |    0.000000000000 ✔
 1 |  2 |  8 |    0.000000000000 |   -0.000000000000 ✔
 1 |  2 |  9 |    0.000000000000 |    0.000000000000 ✔
 1 |  3 |  1 |    0.000000000000 |    0.000000000000 ✔
 1 |  3 |  2 |    0.000000000000 |    0.000000000000 ✔
 1 |  3 |  3 |    3.428571428571 |    3.428571428571 ✔
 1 |  3 |  4 |    0.000000000000 |    0.000000000000 ✔
 1 |  3 |  5 |    0.000000000000 |   -0.000000000000 ✔
 1 |  3 |  6 |    0.000000000000 |   -0.000000000000 ✔
 1 |  3 |  7 |    0.000000000000 |    0.000000000000 ✔
 1 |  3 |  8 |    0.000000000000 |   -0.000000000000 ✔
 1 |  3 |  9 |    0.000000000000 |   -0.000000000000 ✔
 1 |  4 |  1 |    0.000000000000 |    0.000000000000 ✔
 1 |  4 |  2 |    0.000000000000 |    0.000000000000 ✔
 1 |  4 |  3 |    0.000000000000 |    0.000000000000 ✔
 1 |  4 |  4 |    4.444444444444 |    4.444444444444 ✔
 1 |  4 |  5 |    0.000000000000 |    0.000000000000 ✔
 1 |  4 |  6 |    0.000000000000 |    0.000000000000 ✔
 1 |  4 |  7 |    0.000000000000 |    0.000000000000 ✔
 1 |  4 |  8 |    0.000000000000 |   -0.000000000000 ✔
 1 |  4 |  9 |    0.000000000000 |    0.000000000000 ✔
 1 |  5 |  1 |    0.000000000000 |    0.000000000000 ✔
 1 |  5 |  2 |    0.000000000000 |   -0.000000000000 ✔
 1 |  5 |  3 |    0.000000000000 |   -0.000000000000 ✔
 1 |  5 |  4 |    0.000000000000 |    0.000000000000 ✔
 1 |  5 |  5 |    5.454545454545 |    5.454545454545 ✔
 1 |  5 |  6 |    0.000000000000 |    0.000000000000 ✔
 1 |  5 |  7 |    0.000000000000 |   -0.000000000000 ✔
 1 |  5 |  8 |    0.000000000000 |   -0.000000000000 ✔
 1 |  5 |  9 |    0.000000000000 |   -0.000000000000 ✔
 1 |  6 |  1 |    0.000000000000 |   -0.000000000000 ✔
 1 |  6 |  2 |    0.000000000000 |    0.000000000000 ✔
 1 |  6 |  3 |    0.000000000000 |   -0.000000000000 ✔
 1 |  6 |  4 |    0.000000000000 |    0.000000000000 ✔
 1 |  6 |  5 |    0.000000000000 |    0.000000000000 ✔
 1 |  6 |  6 |    6.461538461538 |    6.461538461538 ✔
 1 |  6 |  7 |    0.000000000000 |    0.000000000000 ✔
 1 |  6 |  8 |    0.000000000000 |    0.000000000000 ✔
 1 |  6 |  9 |    0.000000000000 |    0.000000000000 ✔
 1 |  7 |  1 |    0.000000000000 |    0.000000000000 ✔
 1 |  7 |  2 |    0.000000000000 |    0.000000000000 ✔
 1 |  7 |  3 |    0.000000000000 |    0.000000000000 ✔
 1 |  7 |  4 |    0.000000000000 |    0.000000000000 ✔
 1 |  7 |  5 |    0.000000000000 |   -0.000000000000 ✔
 1 |  7 |  6 |    0.000000000000 |    0.000000000000 ✔
 1 |  7 |  7 |    7.466666666667 |    7.466666666667 ✔
 1 |  7 |  8 |    0.000000000000 |    0.000000000000 ✔
 1 |  7 |  9 |    0.000000000000 |    0.000000000000 ✔
 1 |  8 |  1 |    0.000000000000 |    0.000000000000 ✔
 1 |  8 |  2 |    0.000000000000 |   -0.000000000000 ✔
 1 |  8 |  3 |    0.000000000000 |   -0.000000000000 ✔
 1 |  8 |  4 |    0.000000000000 |   -0.000000000000 ✔
 1 |  8 |  5 |    0.000000000000 |   -0.000000000000 ✔
 1 |  8 |  6 |    0.000000000000 |    0.000000000000 ✔
 1 |  8 |  7 |    0.000000000000 |    0.000000000000 ✔
 1 |  8 |  8 |    8.470588235294 |    8.470588235294 ✔
 1 |  8 |  9 |    0.000000000000 |   -0.000000000000 ✔
 1 |  9 |  1 |    0.000000000000 |    0.000000000000 ✔
 1 |  9 |  2 |    0.000000000000 |    0.000000000000 ✔
 1 |  9 |  3 |    0.000000000000 |   -0.000000000000 ✔
 1 |  9 |  4 |    0.000000000000 |    0.000000000000 ✔
 1 |  9 |  5 |    0.000000000000 |   -0.000000000000 ✔
 1 |  9 |  6 |    0.000000000000 |    0.000000000000 ✔
 1 |  9 |  7 |    0.000000000000 |    0.000000000000 ✔
 1 |  9 |  8 |    0.000000000000 |   -0.000000000000 ✔
 1 |  9 |  9 |    9.473684210526 |    9.473684210526 ✔
 2 |  2 |  2 |    9.600000000000 |    9.600000000000 ✔
 2 |  2 |  3 |    0.000000000000 |    0.000000000000 ✔
 2 |  2 |  4 |    0.000000000000 |    0.000000000000 ✔
 2 |  2 |  5 |    0.000000000000 |    0.000000000000 ✔
 2 |  2 |  6 |    0.000000000000 |   -0.000000000000 ✔
 2 |  2 |  7 |    0.000000000000 |    0.000000000000 ✔
 2 |  2 |  8 |    0.000000000000 |    0.000000000000 ✔
 2 |  2 |  9 |    0.000000000000 |    0.000000000000 ✔
 2 |  3 |  2 |    0.000000000000 |    0.000000000000 ✔
 2 |  3 |  3 |   34.285714285714 |   34.285714285714 ✔
 2 |  3 |  4 |    0.000000000000 |    0.000000000000 ✔
 2 |  3 |  5 |    0.000000000000 |    0.000000000000 ✔
 2 |  3 |  6 |    0.000000000000 |    0.000000000000 ✔
 2 |  3 |  7 |    0.000000000000 |   -0.000000000000 ✔
 2 |  3 |  8 |    0.000000000000 |    0.000000000000 ✔
 2 |  3 |  9 |    0.000000000000 |   -0.000000000000 ✔
 2 |  4 |  2 |    0.000000000000 |    0.000000000000 ✔
 2 |  4 |  3 |    0.000000000000 |    0.000000000000 ✔
 2 |  4 |  4 |   80.000000000000 |   80.000000000000 ✔
 2 |  4 |  5 |    0.000000000000 |    0.000000000000 ✔
 2 |  4 |  6 |    0.000000000000 |   -0.000000000000 ✔
 2 |  4 |  7 |    0.000000000000 |   -0.000000000000 ✔
 2 |  4 |  8 |    0.000000000000 |    0.000000000000 ✔
 2 |  4 |  9 |    0.000000000000 |   -0.000000000000 ✔
 2 |  5 |  2 |    0.000000000000 |    0.000000000000 ✔
 2 |  5 |  3 |    0.000000000000 |    0.000000000000 ✔
 2 |  5 |  4 |    0.000000000000 |    0.000000000000 ✔
 2 |  5 |  5 |  152.727272727273 |  152.727272727273 ✔
 2 |  5 |  6 |    0.000000000000 |   -0.000000000000 ✔
 2 |  5 |  7 |    0.000000000000 |   -0.000000000000 ✔
 2 |  5 |  8 |    0.000000000000 |    0.000000000000 ✔
 2 |  5 |  9 |    0.000000000000 |   -0.000000000000 ✔
 2 |  6 |  2 |    0.000000000000 |   -0.000000000000 ✔
 2 |  6 |  3 |    0.000000000000 |    0.000000000000 ✔
 2 |  6 |  4 |    0.000000000000 |   -0.000000000000 ✔
 2 |  6 |  5 |    0.000000000000 |   -0.000000000000 ✔
 2 |  6 |  6 |  258.461538461538 |  258.461538461538 ✔
 2 |  6 |  7 |    0.000000000000 |    0.000000000000 ✔
 2 |  6 |  8 |    0.000000000000 |    0.000000000000 ✔
 2 |  6 |  9 |    0.000000000000 |   -0.000000000000 ✔
 2 |  7 |  2 |    0.000000000000 |    0.000000000000 ✔
 2 |  7 |  3 |    0.000000000000 |   -0.000000000000 ✔
 2 |  7 |  4 |    0.000000000000 |   -0.000000000000 ✔
 2 |  7 |  5 |    0.000000000000 |   -0.000000000000 ✔
 2 |  7 |  6 |    0.000000000000 |    0.000000000000 ✔
 2 |  7 |  7 |  403.200000000000 |  403.200000000000 ✔
 2 |  7 |  8 |    0.000000000000 |   -0.000000000000 ✔
 2 |  7 |  9 |    0.000000000000 |   -0.000000000000 ✔
 2 |  8 |  2 |    0.000000000000 |    0.000000000000 ✔
 2 |  8 |  3 |    0.000000000000 |    0.000000000000 ✔
 2 |  8 |  4 |    0.000000000000 |    0.000000000000 ✔
 2 |  8 |  5 |    0.000000000000 |    0.000000000000 ✔
 2 |  8 |  6 |    0.000000000000 |    0.000000000000 ✔
 2 |  8 |  7 |    0.000000000000 |   -0.000000000000 ✔
 2 |  8 |  8 |  592.941176470588 |  592.941176470588 ✔
 2 |  8 |  9 |    0.000000000000 |    0.000000000000 ✔
 2 |  9 |  2 |    0.000000000000 |    0.000000000000 ✔
 2 |  9 |  3 |    0.000000000000 |   -0.000000000000 ✔
 2 |  9 |  4 |    0.000000000000 |   -0.000000000000 ✔
 2 |  9 |  5 |    0.000000000000 |   -0.000000000000 ✔
 2 |  9 |  6 |    0.000000000000 |   -0.000000000000 ✔
 2 |  9 |  7 |    0.000000000000 |   -0.000000000000 ✔
 2 |  9 |  8 |    0.000000000000 |    0.000000000000 ✔
 2 |  9 |  9 |  833.684210526316 |  833.684210526316 ✔
 3 |  3 |  3 |  205.714285714286 |  205.714285714286 ✔
 3 |  3 |  4 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  5 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  6 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  7 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  8 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  9 |    0.000000000000 |    0.000000000000 ✔
 3 |  4 |  3 |    0.000000000000 |   -0.000000000000 ✔
 3 |  4 |  4 | 1120.000000000000 | 1120.000000000000 ✔
 3 |  4 |  5 |    0.000000000000 |    0.000000000000 ✔
 3 |  4 |  6 |    0.000000000000 |    0.000000000000 ✔
 3 |  4 |  7 |    0.000000000000 |    0.000000000000 ✔
 3 |  4 |  8 |    0.000000000000 |   -0.000000000000 ✔
 3 |  4 |  9 |    0.000000000000 |   -0.000000000000 ✔
 3 |  5 |  3 |    0.000000000000 |   -0.000000000000 ✔
 3 |  5 |  4 |    0.000000000000 |    0.000000000000 ✔
 3 |  5 |  5 | 3665.454545454545 | 3665.454545454545 ✔
 3 |  5 |  6 |    0.000000000000 |    0.000000000000 ✔
 3 |  5 |  7 |    0.000000000000 |   -0.000000000000 ✔
 3 |  5 |  8 |    0.000000000000 |   -0.000000000000 ✔
 3 |  5 |  9 |    0.000000000000 |   -0.000000000000 ✔
 3 |  6 |  3 |    0.000000000000 |    0.000000000000 ✔
 3 |  6 |  4 |    0.000000000000 |    0.000000000000 ✔
 3 |  6 |  5 |    0.000000000000 |    0.000000000000 ✔
 3 |  6 |  6 | 9304.615384615385 | 9304.615384615387 ✔
 3 |  6 |  7 |    0.000000000000 |   -0.000000000000 ✔
 3 |  6 |  8 |    0.000000000000 |   -0.000000000000 ✔
 3 |  6 |  9 |    0.000000000000 |   -0.000000000002 ✔
 3 |  7 |  3 |    0.000000000000 |   -0.000000000000 ✔
 3 |  7 |  4 |    0.000000000000 |    0.000000000000 ✔
 3 |  7 |  5 |    0.000000000000 |   -0.000000000000 ✔
 3 |  7 |  6 |    0.000000000000 |   -0.000000000000 ✔
 3 |  7 |  7 | 20160.000000000000 | 20160.000000000004 ✔
 3 |  7 |  8 |    0.000000000000 |    0.000000000000 ✔
 3 |  7 |  9 |    0.000000000000 |   -0.000000000003 ✔
 3 |  8 |  3 |    0.000000000000 |    0.000000000000 ✔
 3 |  8 |  4 |    0.000000000000 |   -0.000000000000 ✔
 3 |  8 |  5 |    0.000000000000 |   -0.000000000000 ✔
 3 |  8 |  6 |    0.000000000000 |   -0.000000000000 ✔
 3 |  8 |  7 |    0.000000000000 |    0.000000000000 ✔
 3 |  8 |  8 | 39134.117647058825 | 39134.117647058825 ✔
 3 |  8 |  9 |    0.000000000000 |    0.000000000000 ✔
 3 |  9 |  3 |    0.000000000000 |    0.000000000000 ✔
 3 |  9 |  4 |    0.000000000000 |   -0.000000000000 ✔
 3 |  9 |  5 |    0.000000000000 |   -0.000000000000 ✔
 3 |  9 |  6 |    0.000000000000 |   -0.000000000002 ✔
 3 |  9 |  7 |    0.000000000000 |   -0.000000000003 ✔
 3 |  9 |  8 |    0.000000000000 |    0.000000000000 ✔
 3 |  9 |  9 | 70029.473684210534 | 70029.473684210505 ✔
 4 |  4 |  4 | 8960.000000000000 | 8960.000000000002 ✔
 4 |  4 |  5 |    0.000000000000 |   -0.000000000002 ✔
 4 |  4 |  6 |    0.000000000000 |   -0.000000000001 ✔
 4 |  4 |  7 |    0.000000000000 |   -0.000000000000 ✔
 4 |  4 |  8 |    0.000000000000 |    0.000000000007 ✔
 4 |  4 |  9 |    0.000000000000 |    0.000000000000 ✔
 4 |  5 |  4 |    0.000000000000 |   -0.000000000002 ✔
 4 |  5 |  5 | 65978.181818181823 | 65978.181818181838 ✔
 4 |  5 |  6 |    0.000000000000 |   -0.000000000001 ✔
 4 |  5 |  7 |    0.000000000000 |   -0.000000000058 ✔
 4 |  5 |  8 |    0.000000000000 |   -0.000000000002 ✔
 4 |  5 |  9 |    0.000000000000 |   -0.000000000007 ✔
 4 |  6 |  4 |    0.000000000000 |   -0.000000000001 ✔
 4 |  6 |  5 |    0.000000000000 |   -0.000000000001 ✔
 4 |  6 |  6 | 279138.461538461561 | 279138.461538461503 ✔
 4 |  6 |  7 |    0.000000000000 |   -0.000000000018 ✔
 4 |  6 |  8 |    0.000000000000 |    0.000000000055 ✔
 4 |  6 |  9 |    0.000000000000 |    0.000000000029 ✔
 4 |  7 |  4 |    0.000000000000 |   -0.000000000000 ✔
 4 |  7 |  5 |    0.000000000000 |   -0.000000000058 ✔
 4 |  7 |  6 |    0.000000000000 |   -0.000000000018 ✔
 4 |  7 |  7 | 887040.000000000000 | 887040.000000000000 ✔
 4 |  7 |  8 |    0.000000000000 |    0.000000000031 ✔
 4 |  7 |  9 |    0.000000000000 |    0.000000000104 ✔
 4 |  8 |  4 |    0.000000000000 |    0.000000000007 ✔
 4 |  8 |  5 |    0.000000000000 |   -0.000000000002 ✔
 4 |  8 |  6 |    0.000000000000 |    0.000000000055 ✔
 4 |  8 |  7 |    0.000000000000 |    0.000000000031 ✔
 4 |  8 |  8 | 2348047.058823529165 | 2348047.058823529631 ✔
 4 |  8 |  9 |    0.000000000000 |   -0.000000000015 ✔
 4 |  9 |  4 |    0.000000000000 |    0.000000000000 ✔
 4 |  9 |  5 |    0.000000000000 |   -0.000000000007 ✔
 4 |  9 |  6 |    0.000000000000 |    0.000000000029 ✔
 4 |  9 |  7 |    0.000000000000 |    0.000000000104 ✔
 4 |  9 |  8 |    0.000000000000 |   -0.000000000015 ✔
 4 |  9 |  9 | 5462298.947368421592 | 5462298.947368418798 ✔
 5 |  5 |  5 | 659781.818181818235 | 659781.818181818351 ✔
 5 |  5 |  6 |    0.000000000000 |   -0.000000000002 ✔
 5 |  5 |  7 |    0.000000000000 |    0.000000000233 ✔
 5 |  5 |  8 |    0.000000000000 |    0.000000000567 ✔
 5 |  5 |  9 |    0.000000000000 |    0.000000000000 ✔
 5 |  6 |  5 |    0.000000000000 |   -0.000000000002 ✔
 5 |  6 |  6 | 6141046.153846153989 | 6141046.153846156783 ✔
 5 |  6 |  7 |    0.000000000000 |    0.000000000250 ✔
 5 |  6 |  8 |    0.000000000000 |    0.000000001630 ✔
 5 |  6 |  9 |    0.000000000000 |    0.000000000931 ✔
 5 |  7 |  5 |    0.000000000000 |    0.000000000233 ✔
 5 |  7 |  6 |    0.000000000000 |    0.000000000250 ✔
 5 |  7 |  7 | 31933440.000000000000 | 31933440.000000000000 ✔
 5 |  7 |  8 |    0.000000000000 |    0.000000002503 ✔
 5 |  7 |  9 |    0.000000000000 |    0.000000003725 ✔
 5 |  8 |  5 |    0.000000000000 |    0.000000000567 ✔
 5 |  8 |  6 |    0.000000000000 |    0.000000001630 ✔
 5 |  8 |  7 |    0.000000000000 |    0.000000002503 ✔
 5 |  8 |  8 | 122098447.058823525906 | 122098447.058823525906 ✔
 5 |  8 |  9 |    0.000000000000 |   -0.000000001397 ✔
 5 |  9 |  5 |    0.000000000000 |    0.000000000000 ✔
 5 |  9 |  6 |    0.000000000000 |    0.000000000931 ✔
 5 |  9 |  7 |    0.000000000000 |    0.000000003725 ✔
 5 |  9 |  8 |    0.000000000000 |   -0.000000001397 ✔
 5 |  9 |  9 | 382360926.315789461136 | 382360926.315789461136 ✔
Test Summary:                              | Pass  Total  Time
∫Pᵢᵐ(x)Pⱼᵐ(x)dx = 2(j+m)!/(2j+1)(j-m)! δᵢⱼ |  355    355  0.9s
```

#### Normalization & Orthogonality of $Y_{lm}(\theta,\varphi)$

```math
\int_0^{2\pi}
\int_0^\pi
Y_{lm}(\theta,\varphi)^* Y_{l'm'}(\theta,\varphi) \sin(\theta)
~\mathrm{d}\theta \mathrm{d}\varphi
= \delta_{ll'} \delta_{mm'}
```
```
l₁ | l₂ | m₁ | m₂ |        analytical |         numerical 
-- | -- | -- | -- | ----------------- | ----------------- 
 0 |  0 |  0 |  0 |    1.000000000000 |    1.000000000000 ✔
 0 |  1 |  0 | -1 |    0.000000000000 |    0.000000000000 ✔
 0 |  1 |  0 |  0 |    0.000000000000 |   -0.000000000000 ✔
 0 |  1 |  0 |  1 |    0.000000000000 |    0.000000000000 ✔
 0 |  2 |  0 | -2 |    0.000000000000 |   -0.000000000000 ✔
 0 |  2 |  0 | -1 |    0.000000000000 |    0.000000000000 ✔
 0 |  2 |  0 |  0 |    0.000000000000 |    0.000000000000 ✔
 0 |  2 |  0 |  1 |    0.000000000000 |   -0.000000000000 ✔
 0 |  2 |  0 |  2 |    0.000000000000 |   -0.000000000000 ✔
 1 |  0 | -1 |  0 |    0.000000000000 |    0.000000000000 ✔
 1 |  0 |  0 |  0 |    0.000000000000 |   -0.000000000000 ✔
 1 |  0 |  1 |  0 |    0.000000000000 |    0.000000000000 ✔
 1 |  1 | -1 | -1 |    1.000000000000 |    1.000000000000 ✔
 1 |  1 | -1 |  0 |    0.000000000000 |    0.000000000000 ✔
 1 |  1 | -1 |  1 |    0.000000000000 |    0.000000000000 ✔
 1 |  1 |  0 | -1 |    0.000000000000 |    0.000000000000 ✔
 1 |  1 |  0 |  0 |    1.000000000000 |    1.000000000000 ✔
 1 |  1 |  0 |  1 |    0.000000000000 |   -0.000000000000 ✔
 1 |  1 |  1 | -1 |    0.000000000000 |    0.000000000000 ✔
 1 |  1 |  1 |  0 |    0.000000000000 |   -0.000000000000 ✔
 1 |  1 |  1 |  1 |    1.000000000000 |    1.000000000000 ✔
 1 |  2 | -1 | -2 |    0.000000000000 |   -0.000000000000 ✔
 1 |  2 | -1 | -1 |    0.000000000000 |   -0.000000000000 ✔
 1 |  2 | -1 |  0 |    0.000000000000 |    0.000000000000 ✔
 1 |  2 | -1 |  1 |    0.000000000000 |   -0.000000000000 ✔
 1 |  2 | -1 |  2 |    0.000000000000 |    0.000000000000 ✔
 1 |  2 |  0 | -2 |    0.000000000000 |   -0.000000000000 ✔
 1 |  2 |  0 | -1 |    0.000000000000 |   -0.000000000000 ✔
 1 |  2 |  0 |  0 |    0.000000000000 |    0.000000000000 ✔
 1 |  2 |  0 |  1 |    0.000000000000 |    0.000000000000 ✔
 1 |  2 |  0 |  2 |    0.000000000000 |   -0.000000000000 ✔
 1 |  2 |  1 | -2 |    0.000000000000 |   -0.000000000000 ✔
 1 |  2 |  1 | -1 |    0.000000000000 |   -0.000000000000 ✔
 1 |  2 |  1 |  0 |    0.000000000000 |   -0.000000000000 ✔
 1 |  2 |  1 |  1 |    0.000000000000 |   -0.000000000000 ✔
 1 |  2 |  1 |  2 |    0.000000000000 |    0.000000000000 ✔
 2 |  0 | -2 |  0 |    0.000000000000 |   -0.000000000000 ✔
 2 |  0 | -1 |  0 |    0.000000000000 |    0.000000000000 ✔
 2 |  0 |  0 |  0 |    0.000000000000 |    0.000000000000 ✔
 2 |  0 |  1 |  0 |    0.000000000000 |   -0.000000000000 ✔
 2 |  0 |  2 |  0 |    0.000000000000 |   -0.000000000000 ✔
 2 |  1 | -2 | -1 |    0.000000000000 |   -0.000000000000 ✔
 2 |  1 | -2 |  0 |    0.000000000000 |   -0.000000000000 ✔
 2 |  1 | -2 |  1 |    0.000000000000 |   -0.000000000000 ✔
 2 |  1 | -1 | -1 |    0.000000000000 |   -0.000000000000 ✔
 2 |  1 | -1 |  0 |    0.000000000000 |   -0.000000000000 ✔
 2 |  1 | -1 |  1 |    0.000000000000 |   -0.000000000000 ✔
 2 |  1 |  0 | -1 |    0.000000000000 |    0.000000000000 ✔
 2 |  1 |  0 |  0 |    0.000000000000 |    0.000000000000 ✔
 2 |  1 |  0 |  1 |    0.000000000000 |   -0.000000000000 ✔
 2 |  1 |  1 | -1 |    0.000000000000 |   -0.000000000000 ✔
 2 |  1 |  1 |  0 |    0.000000000000 |    0.000000000000 ✔
 2 |  1 |  1 |  1 |    0.000000000000 |   -0.000000000000 ✔
 2 |  1 |  2 | -1 |    0.000000000000 |    0.000000000000 ✔
 2 |  1 |  2 |  0 |    0.000000000000 |   -0.000000000000 ✔
 2 |  1 |  2 |  1 |    0.000000000000 |    0.000000000000 ✔
 2 |  2 | -2 | -2 |    1.000000000000 |    1.000000000000 ✔
 2 |  2 | -2 | -1 |    0.000000000000 |   -0.000000000000 ✔
 2 |  2 | -2 |  0 |    0.000000000000 |    0.000000000000 ✔
 2 |  2 | -2 |  1 |    0.000000000000 |    0.000000000000 ✔
 2 |  2 | -2 |  2 |    0.000000000000 |   -0.000000000000 ✔
 2 |  2 | -1 | -2 |    0.000000000000 |   -0.000000000000 ✔
 2 |  2 | -1 | -1 |    1.000000000000 |    1.000000000000 ✔
 2 |  2 | -1 |  0 |    0.000000000000 |   -0.000000000000 ✔
 2 |  2 | -1 |  1 |    0.000000000000 |    0.000000000000 ✔
 2 |  2 | -1 |  2 |    0.000000000000 |   -0.000000000000 ✔
 2 |  2 |  0 | -2 |    0.000000000000 |    0.000000000000 ✔
 2 |  2 |  0 | -1 |    0.000000000000 |   -0.000000000000 ✔
 2 |  2 |  0 |  0 |    1.000000000000 |    1.000000000000 ✔
 2 |  2 |  0 |  1 |    0.000000000000 |    0.000000000000 ✔
 2 |  2 |  0 |  2 |    0.000000000000 |    0.000000000000 ✔
 2 |  2 |  1 | -2 |    0.000000000000 |    0.000000000000 ✔
 2 |  2 |  1 | -1 |    0.000000000000 |    0.000000000000 ✔
 2 |  2 |  1 |  0 |    0.000000000000 |    0.000000000000 ✔
 2 |  2 |  1 |  1 |    1.000000000000 |    1.000000000000 ✔
 2 |  2 |  1 |  2 |    0.000000000000 |    0.000000000000 ✔
 2 |  2 |  2 | -2 |    0.000000000000 |   -0.000000000000 ✔
 2 |  2 |  2 | -1 |    0.000000000000 |   -0.000000000000 ✔
 2 |  2 |  2 |  0 |    0.000000000000 |    0.000000000000 ✔
 2 |  2 |  2 |  1 |    0.000000000000 |    0.000000000000 ✔
 2 |  2 |  2 |  2 |    1.000000000000 |    1.000000000000 ✔
Test Summary:                              | Pass  Total  Time
∫Yₗ₁ₘ₁(θ,φ)Yₗ₂ₘ₂(θ,φ)sinθdθdφ = δₗ₁ₗ₂δₘ₁ₘ₂ |   81     81  1.8s
```

#### Associated Laguerre Polynomials $L_n^{k}(x)$

```math
  \begin{aligned}
  L_n^{k}(x)
    &= \frac{\mathrm{d}^k}{\mathrm{d}x^k} L_n(x) \\
    &= \frac{\mathrm{d}^k}{\mathrm{d}x^k} \frac{1}{n!} \mathrm{e}^x \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left( \mathrm{e}^{-x} x^n \right) \\
    &= \sum_{m=0}^{n-k} (-1)^{m+k} \frac{n!}{m!(m+k)!(n-m-k)!} x^m \\
    &= (-1)^k L_{n-k}^{(k)}(x)
  \end{aligned}
```

``n=0, k=0:`` ✔
```math
\begin{aligned}
  L_{0}^{0}(x)
   = e^{x} e^{ - x}
  &= 1 \\
  &= 1 \\
  &= 1
\end{aligned}
```

``n=1, k=0:`` ✔
```math
\begin{aligned}
  L_{1}^{0}(x)
   = e^{x} \frac{\mathrm{d}}{\mathrm{d}x} x e^{ - x}
  &= 1 - x \\
  &= 1 - x \\
  &= 1 - x
\end{aligned}
```

``n=1, k=1:`` ✔
```math
\begin{aligned}
  L_{1}^{1}(x)
   = \frac{\mathrm{d}}{\mathrm{d}x} e^{x} \frac{\mathrm{d}}{\mathrm{d}x} x e^{ - x}
  &= -1 \\
  &= -1 \\
  &= -1
\end{aligned}
```

``n=2, k=0:`` ✔
```math
\begin{aligned}
  L_{2}^{0}(x)
   = \frac{1}{2} e^{x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} x^{2} e^{ - x}
  &= 1 - 2 x + \frac{1}{2} x^{2} \\
  &= 1 - 2 x + \frac{1}{2} x^{2} \\
  &= 1 - 2 x + \frac{1}{2} x^{2}
\end{aligned}
```

``n=2, k=1:`` ✔
```math
\begin{aligned}
  L_{2}^{1}(x)
   = \frac{\mathrm{d}}{\mathrm{d}x} \frac{1}{2} e^{x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} x^{2} e^{ - x}
  &= -2 + x \\
  &= -2 + x \\
  &= -2 + x
\end{aligned}
```

``n=2, k=2:`` ✔
```math
\begin{aligned}
  L_{2}^{2}(x)
   = \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{1}{2} e^{x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} x^{2} e^{ - x}
  &= 1 \\
  &= 1 \\
  &= 1
\end{aligned}
```

``n=3, k=0:`` ✔
```math
\begin{aligned}
  L_{3}^{0}(x)
   = \frac{1}{6} e^{x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} x^{3} e^{ - x}
  &= 1 - \frac{1}{6} x^{3} - 3 x + \frac{3}{2} x^{2} \\
  &= 1 - \frac{1}{6} x^{3} - 3 x + \frac{3}{2} x^{2} \\
  &= 1 - \frac{1}{6} x^{3} - 3 x + \frac{3}{2} x^{2}
\end{aligned}
```

``n=3, k=1:`` ✔
```math
\begin{aligned}
  L_{3}^{1}(x)
   = \frac{\mathrm{d}}{\mathrm{d}x} \frac{1}{6} e^{x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} x^{3} e^{ - x}
  &= -3 + 3 x - \frac{1}{2} x^{2} \\
  &= -3 + 3 x - \frac{1}{2} x^{2} \\
  &= -3 + 3 x - \frac{1}{2} x^{2}
\end{aligned}
```

``n=3, k=2:`` ✔
```math
\begin{aligned}
  L_{3}^{2}(x)
   = \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{1}{6} e^{x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} x^{3} e^{ - x}
  &= 3 - x \\
  &= 3 - x \\
  &= 3 - x
\end{aligned}
```

``n=3, k=3:`` ✔
```math
\begin{aligned}
  L_{3}^{3}(x)
   = \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{1}{6} e^{x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} x^{3} e^{ - x}
  &= -1 \\
  &= -1 \\
  &= -1
\end{aligned}
```

``n=4, k=0:`` ✔
```math
\begin{aligned}
  L_{4}^{0}(x)
   = \frac{1}{24} e^{x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} x^{4} e^{ - x}
  &= 1 - \frac{2}{3} x^{3} - 4 x + 3 x^{2} + \frac{1}{24} x^{4} \\
  &= 1 - \frac{2}{3} x^{3} - 4 x + 3 x^{2} + \frac{1}{24} x^{4} \\
  &= 1 - \frac{2}{3} x^{3} - 4 x + 3 x^{2} + \frac{1}{24} x^{4}
\end{aligned}
```

``n=4, k=1:`` ✔
```math
\begin{aligned}
  L_{4}^{1}(x)
   = \frac{\mathrm{d}}{\mathrm{d}x} \frac{1}{24} e^{x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} x^{4} e^{ - x}
  &= -4 + \frac{1}{6} x^{3} + 6 x - 2 x^{2} \\
  &= -4 + \frac{1}{6} x^{3} + 6 x - 2 x^{2} \\
  &= -4 + \frac{1}{6} x^{3} + 6 x - 2 x^{2}
\end{aligned}
```

``n=4, k=2:`` ✔
```math
\begin{aligned}
  L_{4}^{2}(x)
   = \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{1}{24} e^{x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} x^{4} e^{ - x}
  &= 6 - 4 x + \frac{1}{2} x^{2} \\
  &= 6 - 4 x + \frac{1}{2} x^{2} \\
  &= 6 - 4 x + \frac{1}{2} x^{2}
\end{aligned}
```

``n=4, k=3:`` ✔
```math
\begin{aligned}
  L_{4}^{3}(x)
   = \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{1}{24} e^{x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} x^{4} e^{ - x}
  &= -4 + x \\
  &= -4 + x \\
  &= -4 + x
\end{aligned}
```

``n=4, k=4:`` ✔
```math
\begin{aligned}
  L_{4}^{4}(x)
   = \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{1}{24} e^{x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x} x^{4} e^{ - x}
  &= 1 \\
  &= 1 \\
  &= 1
\end{aligned}
```

```
Test Summary:                                          | Pass  Total  Time
Lₙᵏ(x) = dᵏ/dxᵏ Lₙ(x); Lₙ(x) = 1/(n!) eˣ dⁿ/dxⁿ e⁻ˣ xⁿ |   15     15  0.5s
```

#### Normalization & Orthogonality of $L_n^{k}(x)$

```math
\int_{0}^{\infty} \mathrm{e}^{-x} x^k L_i^k(x) L_j^k(x) \mathrm{d}x = \frac{i!}{(i-k)!} \delta_{ij}
```

Replace $n+k$ with $n$ for [the definition of Wolfram MathWorld](https://mathworld.wolfram.com/AssociatedLaguerrePolynomial.html).
```
 i |  j |  k |        analytical |         numerical 
-- | -- | -- | ----------------- | ----------------- 
 0 |  0 |  0 |    1.000000000000 |    1.000000000000 ✔
 0 |  1 |  0 |    0.000000000000 |    0.000000000000 ✔
 0 |  2 |  0 |    0.000000000000 |    0.000000000000 ✔
 0 |  3 |  0 |    0.000000000000 |    0.000000000000 ✔
 0 |  4 |  0 |    0.000000000000 |    0.000000000000 ✔
 0 |  5 |  0 |    0.000000000000 |   -0.000000000000 ✔
 0 |  6 |  0 |    0.000000000000 |   -0.000000000000 ✔
 0 |  7 |  0 |    0.000000000000 |    0.000000000000 ✔
 1 |  0 |  0 |    0.000000000000 |    0.000000000000 ✔
 1 |  1 |  0 |    1.000000000000 |    1.000000000000 ✔
 1 |  1 |  1 |    1.000000000000 |    1.000000000000 ✔
 1 |  2 |  0 |    0.000000000000 |   -0.000000000000 ✔
 1 |  2 |  1 |    0.000000000000 |    0.000000000000 ✔
 1 |  3 |  0 |    0.000000000000 |   -0.000000000000 ✔
 1 |  3 |  1 |    0.000000000000 |   -0.000000000000 ✔
 1 |  4 |  0 |    0.000000000000 |   -0.000000000000 ✔
 1 |  4 |  1 |    0.000000000000 |    0.000000000000 ✔
 1 |  5 |  0 |    0.000000000000 |    0.000000000000 ✔
 1 |  5 |  1 |    0.000000000000 |   -0.000000000000 ✔
 1 |  6 |  0 |    0.000000000000 |    0.000000000000 ✔
 1 |  6 |  1 |    0.000000000000 |   -0.000000000000 ✔
 1 |  7 |  0 |    0.000000000000 |   -0.000000000000 ✔
 1 |  7 |  1 |    0.000000000000 |    0.000000000000 ✔
 2 |  0 |  0 |    0.000000000000 |    0.000000000000 ✔
 2 |  1 |  0 |    0.000000000000 |    0.000000000000 ✔
 2 |  1 |  1 |    0.000000000000 |    0.000000000000 ✔
 2 |  2 |  0 |    1.000000000000 |    1.000000000000 ✔
 2 |  2 |  1 |    2.000000000000 |    2.000000000000 ✔
 2 |  2 |  2 |    2.000000000000 |    2.000000000000 ✔
 2 |  3 |  0 |    0.000000000000 |    0.000000000000 ✔
 2 |  3 |  1 |    0.000000000000 |   -0.000000000000 ✔
 2 |  3 |  2 |    0.000000000000 |   -0.000000000000 ✔
 2 |  4 |  0 |    0.000000000000 |    0.000000000000 ✔
 2 |  4 |  1 |    0.000000000000 |   -0.000000000000 ✔
 2 |  4 |  2 |    0.000000000000 |   -0.000000000000 ✔
 2 |  5 |  0 |    0.000000000000 |    0.000000000000 ✔
 2 |  5 |  1 |    0.000000000000 |    0.000000000000 ✔
 2 |  5 |  2 |    0.000000000000 |    0.000000000000 ✔
 2 |  6 |  0 |    0.000000000000 |   -0.000000000000 ✔
 2 |  6 |  1 |    0.000000000000 |    0.000000000000 ✔
 2 |  6 |  2 |    0.000000000000 |   -0.000000000000 ✔
 2 |  7 |  0 |    0.000000000000 |    0.000000000000 ✔
 2 |  7 |  1 |    0.000000000000 |   -0.000000000000 ✔
 2 |  7 |  2 |    0.000000000000 |    0.000000000000 ✔
 3 |  0 |  0 |    0.000000000000 |    0.000000000000 ✔
 3 |  1 |  0 |    0.000000000000 |   -0.000000000000 ✔
 3 |  1 |  1 |    0.000000000000 |   -0.000000000000 ✔
 3 |  2 |  0 |    0.000000000000 |    0.000000000000 ✔
 3 |  2 |  1 |    0.000000000000 |   -0.000000000000 ✔
 3 |  2 |  2 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  0 |    1.000000000000 |    1.000000000000 ✔
 3 |  3 |  1 |    3.000000000000 |    3.000000000000 ✔
 3 |  3 |  2 |    6.000000000000 |    6.000000000000 ✔
 3 |  3 |  3 |    6.000000000000 |    6.000000000000 ✔
 3 |  4 |  0 |    0.000000000000 |    0.000000000000 ✔
 3 |  4 |  1 |    0.000000000000 |    0.000000000000 ✔
 3 |  4 |  2 |    0.000000000000 |   -0.000000000000 ✔
 3 |  4 |  3 |    0.000000000000 |   -0.000000000000 ✔
 3 |  5 |  0 |    0.000000000000 |   -0.000000000000 ✔
 3 |  5 |  1 |    0.000000000000 |   -0.000000000000 ✔
 3 |  5 |  2 |    0.000000000000 |   -0.000000000000 ✔
 3 |  5 |  3 |    0.000000000000 |    0.000000000000 ✔
 3 |  6 |  0 |    0.000000000000 |    0.000000000000 ✔
 3 |  6 |  1 |    0.000000000000 |   -0.000000000000 ✔
 3 |  6 |  2 |    0.000000000000 |    0.000000000000 ✔
 3 |  6 |  3 |    0.000000000000 |    0.000000000000 ✔
 3 |  7 |  0 |    0.000000000000 |   -0.000000000000 ✔
 3 |  7 |  1 |    0.000000000000 |    0.000000000000 ✔
 3 |  7 |  2 |    0.000000000000 |   -0.000000000000 ✔
 3 |  7 |  3 |    0.000000000000 |   -0.000000000000 ✔
 4 |  0 |  0 |    0.000000000000 |    0.000000000000 ✔
 4 |  1 |  0 |    0.000000000000 |   -0.000000000000 ✔
 4 |  1 |  1 |    0.000000000000 |    0.000000000000 ✔
 4 |  2 |  0 |    0.000000000000 |    0.000000000000 ✔
 4 |  2 |  1 |    0.000000000000 |   -0.000000000000 ✔
 4 |  2 |  2 |    0.000000000000 |   -0.000000000000 ✔
 4 |  3 |  0 |    0.000000000000 |    0.000000000000 ✔
 4 |  3 |  1 |    0.000000000000 |    0.000000000000 ✔
 4 |  3 |  2 |    0.000000000000 |    0.000000000000 ✔
 4 |  3 |  3 |    0.000000000000 |   -0.000000000000 ✔
 4 |  4 |  0 |    1.000000000000 |    1.000000000000 ✔
 4 |  4 |  1 |    4.000000000000 |    4.000000000000 ✔
 4 |  4 |  2 |   12.000000000000 |   12.000000000000 ✔
 4 |  4 |  3 |   24.000000000000 |   24.000000000000 ✔
 4 |  4 |  4 |   24.000000000000 |   24.000000000000 ✔
 4 |  5 |  0 |    0.000000000000 |    0.000000000000 ✔
 4 |  5 |  1 |    0.000000000000 |    0.000000000000 ✔
 4 |  5 |  2 |    0.000000000000 |    0.000000000000 ✔
 4 |  5 |  3 |    0.000000000000 |    0.000000000000 ✔
 4 |  5 |  4 |    0.000000000000 |   -0.000000000000 ✔
 4 |  6 |  0 |    0.000000000000 |   -0.000000000000 ✔
 4 |  6 |  1 |    0.000000000000 |    0.000000000000 ✔
 4 |  6 |  2 |    0.000000000000 |   -0.000000000000 ✔
 4 |  6 |  3 |    0.000000000000 |   -0.000000000000 ✔
 4 |  6 |  4 |    0.000000000000 |    0.000000000000 ✔
 4 |  7 |  0 |    0.000000000000 |    0.000000000000 ✔
 4 |  7 |  1 |    0.000000000000 |   -0.000000000000 ✔
 4 |  7 |  2 |    0.000000000000 |    0.000000000000 ✔
 4 |  7 |  3 |    0.000000000000 |    0.000000000000 ✔
 4 |  7 |  4 |    0.000000000000 |    0.000000000000 ✔
 5 |  0 |  0 |    0.000000000000 |   -0.000000000000 ✔
 5 |  1 |  0 |    0.000000000000 |    0.000000000000 ✔
 5 |  1 |  1 |    0.000000000000 |   -0.000000000000 ✔
 5 |  2 |  0 |    0.000000000000 |    0.000000000000 ✔
 5 |  2 |  1 |    0.000000000000 |    0.000000000000 ✔
 5 |  2 |  2 |    0.000000000000 |    0.000000000000 ✔
 5 |  3 |  0 |    0.000000000000 |   -0.000000000000 ✔
 5 |  3 |  1 |    0.000000000000 |   -0.000000000000 ✔
 5 |  3 |  2 |    0.000000000000 |   -0.000000000000 ✔
 5 |  3 |  3 |    0.000000000000 |    0.000000000000 ✔
 5 |  4 |  0 |    0.000000000000 |    0.000000000000 ✔
 5 |  4 |  1 |    0.000000000000 |    0.000000000000 ✔
 5 |  4 |  2 |    0.000000000000 |    0.000000000000 ✔
 5 |  4 |  3 |    0.000000000000 |    0.000000000000 ✔
 5 |  4 |  4 |    0.000000000000 |   -0.000000000000 ✔
 5 |  5 |  0 |    1.000000000000 |    1.000000000000 ✔
 5 |  5 |  1 |    5.000000000000 |    4.999999999999 ✔
 5 |  5 |  2 |   20.000000000000 |   20.000000000000 ✔
 5 |  5 |  3 |   60.000000000000 |   60.000000000000 ✔
 5 |  5 |  4 |  120.000000000000 |  120.000000000000 ✔
 5 |  5 |  5 |  120.000000000000 |  120.000000000000 ✔
 5 |  6 |  0 |    0.000000000000 |    0.000000000000 ✔
 5 |  6 |  1 |    0.000000000000 |   -0.000000000000 ✔
 5 |  6 |  2 |    0.000000000000 |    0.000000000000 ✔
 5 |  6 |  3 |    0.000000000000 |    0.000000000000 ✔
 5 |  6 |  4 |    0.000000000000 |    0.000000000000 ✔
 5 |  6 |  5 |    0.000000000000 |    0.000000000000 ✔
 5 |  7 |  0 |    0.000000000000 |   -0.000000000000 ✔
 5 |  7 |  1 |    0.000000000000 |   -0.000000000000 ✔
 5 |  7 |  2 |    0.000000000000 |   -0.000000000000 ✔
 5 |  7 |  3 |    0.000000000000 |   -0.000000000000 ✔
 5 |  7 |  4 |    0.000000000000 |   -0.000000000000 ✔
 5 |  7 |  5 |    0.000000000000 |   -0.000000000000 ✔
 6 |  0 |  0 |    0.000000000000 |   -0.000000000000 ✔
 6 |  1 |  0 |    0.000000000000 |    0.000000000000 ✔
 6 |  1 |  1 |    0.000000000000 |   -0.000000000000 ✔
 6 |  2 |  0 |    0.000000000000 |   -0.000000000000 ✔
 6 |  2 |  1 |    0.000000000000 |    0.000000000000 ✔
 6 |  2 |  2 |    0.000000000000 |   -0.000000000000 ✔
 6 |  3 |  0 |    0.000000000000 |    0.000000000000 ✔
 6 |  3 |  1 |    0.000000000000 |   -0.000000000000 ✔
 6 |  3 |  2 |    0.000000000000 |    0.000000000000 ✔
 6 |  3 |  3 |    0.000000000000 |    0.000000000000 ✔
 6 |  4 |  0 |    0.000000000000 |   -0.000000000000 ✔
 6 |  4 |  1 |    0.000000000000 |    0.000000000000 ✔
 6 |  4 |  2 |    0.000000000000 |   -0.000000000000 ✔
 6 |  4 |  3 |    0.000000000000 |   -0.000000000000 ✔
 6 |  4 |  4 |    0.000000000000 |    0.000000000000 ✔
 6 |  5 |  0 |    0.000000000000 |    0.000000000000 ✔
 6 |  5 |  1 |    0.000000000000 |   -0.000000000000 ✔
 6 |  5 |  2 |    0.000000000000 |    0.000000000000 ✔
 6 |  5 |  3 |    0.000000000000 |    0.000000000000 ✔
 6 |  5 |  4 |    0.000000000000 |    0.000000000000 ✔
 6 |  5 |  5 |    0.000000000000 |    0.000000000000 ✔
 6 |  6 |  0 |    1.000000000000 |    1.000000000000 ✔
 6 |  6 |  1 |    6.000000000000 |    6.000000000000 ✔
 6 |  6 |  2 |   30.000000000000 |   30.000000000000 ✔
 6 |  6 |  3 |  120.000000000000 |  119.999999999978 ✔
 6 |  6 |  4 |  360.000000000000 |  359.999999999996 ✔
 6 |  6 |  5 |  720.000000000000 |  720.000000000000 ✔
 6 |  6 |  6 |  720.000000000000 |  720.000000000000 ✔
 6 |  7 |  0 |    0.000000000000 |    0.000000000000 ✔
 6 |  7 |  1 |    0.000000000000 |    0.000000000000 ✔
 6 |  7 |  2 |    0.000000000000 |   -0.000000000000 ✔
 6 |  7 |  3 |    0.000000000000 |    0.000000000000 ✔
 6 |  7 |  4 |    0.000000000000 |    0.000000000000 ✔
 6 |  7 |  5 |    0.000000000000 |    0.000000000000 ✔
 6 |  7 |  6 |    0.000000000000 |    0.000000000000 ✔
 7 |  0 |  0 |    0.000000000000 |    0.000000000000 ✔
 7 |  1 |  0 |    0.000000000000 |   -0.000000000000 ✔
 7 |  1 |  1 |    0.000000000000 |    0.000000000000 ✔
 7 |  2 |  0 |    0.000000000000 |    0.000000000000 ✔
 7 |  2 |  1 |    0.000000000000 |   -0.000000000000 ✔
 7 |  2 |  2 |    0.000000000000 |    0.000000000000 ✔
 7 |  3 |  0 |    0.000000000000 |   -0.000000000000 ✔
 7 |  3 |  1 |    0.000000000000 |    0.000000000000 ✔
 7 |  3 |  2 |    0.000000000000 |   -0.000000000000 ✔
 7 |  3 |  3 |    0.000000000000 |   -0.000000000000 ✔
 7 |  4 |  0 |    0.000000000000 |    0.000000000000 ✔
 7 |  4 |  1 |    0.000000000000 |   -0.000000000000 ✔
 7 |  4 |  2 |    0.000000000000 |    0.000000000000 ✔
 7 |  4 |  3 |    0.000000000000 |    0.000000000000 ✔
 7 |  4 |  4 |    0.000000000000 |    0.000000000000 ✔
 7 |  5 |  0 |    0.000000000000 |   -0.000000000000 ✔
 7 |  5 |  1 |    0.000000000000 |   -0.000000000000 ✔
 7 |  5 |  2 |    0.000000000000 |   -0.000000000000 ✔
 7 |  5 |  3 |    0.000000000000 |   -0.000000000000 ✔
 7 |  5 |  4 |    0.000000000000 |   -0.000000000000 ✔
 7 |  5 |  5 |    0.000000000000 |   -0.000000000000 ✔
 7 |  6 |  0 |    0.000000000000 |    0.000000000000 ✔
 7 |  6 |  1 |    0.000000000000 |   -0.000000000000 ✔
 7 |  6 |  2 |    0.000000000000 |   -0.000000000000 ✔
 7 |  6 |  3 |    0.000000000000 |    0.000000000000 ✔
 7 |  6 |  4 |    0.000000000000 |    0.000000000000 ✔
 7 |  6 |  5 |    0.000000000000 |   -0.000000000000 ✔
 7 |  6 |  6 |    0.000000000000 |    0.000000000000 ✔
 7 |  7 |  0 |    1.000000000000 |    1.000000000000 ✔
 7 |  7 |  1 |    7.000000000000 |    7.000000000000 ✔
 7 |  7 |  2 |   42.000000000000 |   42.000000000000 ✔
 7 |  7 |  3 |  210.000000000000 |  210.000000000000 ✔
 7 |  7 |  4 |  840.000000000000 |  840.000000000000 ✔
 7 |  7 |  5 | 2520.000000000000 | 2519.999999999775 ✔
 7 |  7 |  6 | 5040.000000000000 | 5039.999999999985 ✔
 7 |  7 |  7 | 5040.000000000000 | 5040.000000000000 ✔
Test Summary:                                 | Pass  Total  Time
∫exp(-x)xᵏLᵢᵏ(x)Lⱼᵏ(x)dx = (2i+k)!/(i+k)! δᵢⱼ |  204    204  0.8s
```

#### Normalization of $R_{nl}(r)$

```math
\int |R_{nl}(r)|^2 r^2 \mathrm{d}r = 1
```
```
 n |  l |        analytical |         numerical 
-- | -- | ----------------- | ----------------- 
 1 |  0 |    1.000000000000 |    1.000000000000 ✔
 2 |  0 |    1.000000000000 |    1.000000000000 ✔
 2 |  1 |    1.000000000000 |    1.000000000000 ✔
 3 |  0 |    1.000000000000 |    1.000000000000 ✔
 3 |  1 |    1.000000000000 |    0.999999999999 ✔
 3 |  2 |    1.000000000000 |    1.000000000000 ✔
 4 |  0 |    1.000000000000 |    1.000000000000 ✔
 4 |  1 |    1.000000000000 |    1.000000000000 ✔
 4 |  2 |    1.000000000000 |    1.000000000000 ✔
 4 |  3 |    1.000000000000 |    1.000000000000 ✔
 5 |  0 |    1.000000000000 |    1.000000000000 ✔
 5 |  1 |    1.000000000000 |    1.000000000000 ✔
 5 |  2 |    1.000000000000 |    1.000000000000 ✔
 5 |  3 |    1.000000000000 |    1.000000000000 ✔
 5 |  4 |    1.000000000000 |    1.000000000000 ✔
 6 |  0 |    1.000000000000 |    1.000000000000 ✔
 6 |  1 |    1.000000000000 |    1.000000000000 ✔
 6 |  2 |    1.000000000000 |    1.000000000000 ✔
 6 |  3 |    1.000000000000 |    1.000000000000 ✔
 6 |  4 |    1.000000000000 |    1.000000000000 ✔
 6 |  5 |    1.000000000000 |    1.000000000000 ✔
 7 |  0 |    1.000000000000 |    1.000000000000 ✔
 7 |  1 |    1.000000000000 |    1.000000000000 ✔
 7 |  2 |    1.000000000000 |    1.000000000000 ✔
 7 |  3 |    1.000000000000 |    1.000000000000 ✔
 7 |  4 |    1.000000000000 |    1.000000000000 ✔
 7 |  5 |    1.000000000000 |    1.000000000000 ✔
 7 |  6 |    1.000000000000 |    1.000000000000 ✔
 8 |  0 |    1.000000000000 |    1.000000000000 ✔
 8 |  1 |    1.000000000000 |    1.000000000000 ✔
 8 |  2 |    1.000000000000 |    1.000000000000 ✔
 8 |  3 |    1.000000000000 |    1.000000000000 ✔
 8 |  4 |    1.000000000000 |    1.000000000000 ✔
 8 |  5 |    1.000000000000 |    1.000000000000 ✔
 8 |  6 |    1.000000000000 |    1.000000000000 ✔
 8 |  7 |    1.000000000000 |    1.000000000000 ✔
 9 |  0 |    1.000000000000 |    1.000000000000 ✔
 9 |  1 |    1.000000000000 |    1.000000000000 ✔
 9 |  2 |    1.000000000000 |    1.000000000000 ✔
 9 |  3 |    1.000000000000 |    1.000000000000 ✔
 9 |  4 |    1.000000000000 |    1.000000000000 ✔
 9 |  5 |    1.000000000000 |    1.000000000000 ✔
 9 |  6 |    1.000000000000 |    1.000000000000 ✔
 9 |  7 |    1.000000000000 |    1.000000000000 ✔
 9 |  8 |    1.000000000000 |    1.000000000000 ✔
Test Summary:               | Pass  Total  Time
∫|Rₙₗ(r)|²r²dr = δₙ₁ₙ₂δₗ₁ₗ₂ |   45     45  0.6s
```

#### Expected Value of $r$

```math
\langle r \rangle
= \int r |R_{n_1 l_1}(r)|^2 r^2 \mathrm{d}r
= \frac{a_\mu}{2Z} \left[ 3n^2 - l(l+1) \right] \\
a_\mu = a_0 \frac{m_\mathrm{e}}{\mu} \\
\frac{1}{\mu} = \frac{1}{m_\mathrm{e}} + \frac{1}{m_\mathrm{p}}
```

Reference:
- [高柳和夫『朝倉物理学大系 11 原子分子物理学』(2000, 朝倉書店) pp.11-22](https://www.asakura.co.jp/detail.php?book_code=13681)
- [ Quan­tum Me­chan­ics for En­gi­neers by Leon van Dom­me­len](https://web1.eng.famu.fsu.edu/~dommelen/quantum/style_a/nt_rsexp.html)

```
 n |  l |        analytical |         numerical 
-- | -- | ----------------- | ----------------- 
 1 |  0 |    1.500000000000 |    1.500000000000 ✔
 2 |  0 |    6.000000000000 |    6.000000000000 ✔
 2 |  1 |    5.000000000000 |    5.000000000000 ✔
 3 |  0 |   13.500000000000 |   13.500000000000 ✔
 3 |  1 |   12.500000000000 |   12.500000000000 ✔
 3 |  2 |   10.500000000000 |   10.500000000000 ✔
 4 |  0 |   24.000000000000 |   23.999999999999 ✔
 4 |  1 |   23.000000000000 |   22.999999999999 ✔
 4 |  2 |   21.000000000000 |   21.000000000000 ✔
 4 |  3 |   18.000000000000 |   18.000000000000 ✔
 5 |  0 |   37.500000000000 |   37.500000000000 ✔
 5 |  1 |   36.500000000000 |   36.500000000000 ✔
 5 |  2 |   34.500000000000 |   34.500000000000 ✔
 5 |  3 |   31.500000000000 |   31.500000000000 ✔
 5 |  4 |   27.500000000000 |   27.499999999943 ✔
 6 |  0 |   54.000000000000 |   54.000000000001 ✔
 6 |  1 |   53.000000000000 |   53.000000000001 ✔
 6 |  2 |   51.000000000000 |   51.000000000000 ✔
 6 |  3 |   48.000000000000 |   48.000000000000 ✔
 6 |  4 |   44.000000000000 |   44.000000000000 ✔
 6 |  5 |   39.000000000000 |   39.000000000000 ✔
 7 |  0 |   73.500000000000 |   73.500000000000 ✔
 7 |  1 |   72.500000000000 |   72.500000000000 ✔
 7 |  2 |   70.500000000000 |   70.500000000000 ✔
 7 |  3 |   67.500000000000 |   67.500000000000 ✔
 7 |  4 |   63.500000000000 |   63.500000000000 ✔
 7 |  5 |   58.500000000000 |   58.500000000000 ✔
 7 |  6 |   52.500000000000 |   52.499999999992 ✔
 8 |  0 |   96.000000000000 |   96.000000000001 ✔
 8 |  1 |   95.000000000000 |   94.999999999999 ✔
 8 |  2 |   93.000000000000 |   93.000000000000 ✔
 8 |  3 |   90.000000000000 |   90.000000000000 ✔
 8 |  4 |   86.000000000000 |   86.000000000000 ✔
 8 |  5 |   81.000000000000 |   81.000000000000 ✔
 8 |  6 |   75.000000000000 |   75.000000000000 ✔
 8 |  7 |   68.000000000000 |   68.000000000000 ✔
 9 |  0 |  121.500000000000 |  121.500000000001 ✔
 9 |  1 |  120.500000000000 |  120.500000000000 ✔
 9 |  2 |  118.500000000000 |  118.500000000001 ✔
 9 |  3 |  115.500000000000 |  115.500000000000 ✔
 9 |  4 |  111.500000000000 |  111.499999999998 ✔
 9 |  5 |  106.500000000000 |  106.499999999999 ✔
 9 |  6 |  100.500000000000 |  100.500000000000 ✔
 9 |  7 |   93.500000000000 |   93.500000000000 ✔
 9 |  8 |   85.500000000000 |   85.500000000000 ✔
Test Summary:                                                    | Pass  Total  Time
∫r|Rₙₗ(r)|²r²dr = (a₀×mₑ/μ)/2Z × [3n²-l(l+1)]; 1/μ = 1/mₑ + 1/mₚ |   45     45  0.6s
```

#### Expected Value of $r^2$

```math
\langle r^2 \rangle
= \int r^2 |R_{n_1 l_1}(r)|^2 r^2 \mathrm{d}r
= \frac{a_\mu^2}{2Z^2} n^2 \left[ 5n^2 + 1 - 3l(l+1) \right] \\
a_\mu = a_0 \frac{m_\mathrm{e}}{\mu} \\
\frac{1}{\mu} = \frac{1}{m_\mathrm{e}} + \frac{1}{m_\mathrm{p}}
```

Reference:
- [高柳和夫『朝倉物理学大系 11 原子分子物理学』(2000, 朝倉書店) pp.11-22](https://www.asakura.co.jp/detail.php?book_code=13681)
- [ Quan­tum Me­chan­ics for En­gi­neers by Leon van Dom­me­len](https://web1.eng.famu.fsu.edu/~dommelen/quantum/style_a/nt_rsexp.html)
```
 n |  l |        analytical |         numerical 
-- | -- | ----------------- | ----------------- 
 1 |  0 |    3.000000000000 |    3.000000000000 ✔
 2 |  0 |   42.000000000000 |   42.000000000000 ✔
 2 |  1 |   30.000000000000 |   30.000000000000 ✔
 3 |  0 |  207.000000000000 |  207.000000000000 ✔
 3 |  1 |  180.000000000000 |  180.000000000000 ✔
 3 |  2 |  126.000000000000 |  126.000000000000 ✔
 4 |  0 |  648.000000000000 |  647.999999999903 ✔
 4 |  1 |  600.000000000000 |  599.999999999936 ✔
 4 |  2 |  504.000000000000 |  503.999999999975 ✔
 4 |  3 |  360.000000000000 |  359.999999999996 ✔
 5 |  0 | 1575.000000000000 | 1574.999999999999 ✔
 5 |  1 | 1500.000000000000 | 1499.999999999998 ✔
 5 |  2 | 1350.000000000000 | 1350.000000000000 ✔
 5 |  3 | 1125.000000000000 | 1125.000000000003 ✔
 5 |  4 |  825.000000000000 |  825.000000000000 ✔
 6 |  0 | 3258.000000000000 | 3257.999999999997 ✔
 6 |  1 | 3150.000000000000 | 3149.999999999992 ✔
 6 |  2 | 2934.000000000000 | 2933.999999999998 ✔
 6 |  3 | 2610.000000000000 | 2610.000000000033 ✔
 6 |  4 | 2178.000000000000 | 2178.000000000008 ✔
 6 |  5 | 1638.000000000000 | 1638.000000000000 ✔
 7 |  0 | 6027.000000000000 | 6026.999999999992 ✔
 7 |  1 | 5880.000000000000 | 5880.000000000003 ✔
 7 |  2 | 5586.000000000000 | 5585.999999999990 ✔
 7 |  3 | 5145.000000000000 | 5144.999999999992 ✔
 7 |  4 | 4557.000000000000 | 4556.999999999997 ✔
 7 |  5 | 3822.000000000000 | 3821.999999999999 ✔
 7 |  6 | 2940.000000000000 | 2940.000000000001 ✔
 8 |  0 | 10272.000000000000 | 10272.000000000029 ✔
 8 |  1 | 10080.000000000000 | 10079.999999999995 ✔
 8 |  2 | 9696.000000000000 | 9695.999999999993 ✔
 8 |  3 | 9120.000000000000 | 9120.000000000011 ✔
 8 |  4 | 8352.000000000000 | 8352.000000000002 ✔
 8 |  5 | 7392.000000000000 | 7392.000000000010 ✔
 8 |  6 | 6240.000000000000 | 6240.000000000000 ✔
 8 |  7 | 4896.000000000000 | 4896.000000000008 ✔
 9 |  0 | 16443.000000000000 | 16443.000000000102 ✔
 9 |  1 | 16200.000000000000 | 16200.000000000040 ✔
 9 |  2 | 15714.000000000000 | 15714.000000000149 ✔
 9 |  3 | 14985.000000000000 | 14984.999999999918 ✔
 9 |  4 | 14013.000000000000 | 14012.999999999545 ✔
 9 |  5 | 12798.000000000000 | 12797.999999999807 ✔
 9 |  6 | 11340.000000000000 | 11339.999999999945 ✔
 9 |  7 | 9639.000000000000 | 9638.999999999991 ✔
 9 |  8 | 7695.000000000000 | 7694.999999999998 ✔
Test Summary:                                                            | Pass  Total  Time
∫r²|Rₙₗ(r)|²r²dr = (a₀×mₑ/μ)²/2Z² × n²[5n²+1-3l(l+1)]; 1/μ = 1/mₑ + 1/mₚ |   45     45  0.6s
```

#### Virial Theorem

The virial theorem $2\langle T \rangle + \langle V \rangle = 0$ and the definition of Hamiltonian $\langle H \rangle = \langle T \rangle + \langle V \rangle$ derive $\langle H \rangle = \frac{1}{2} \langle V \rangle$ and $\langle H \rangle = -\langle T \rangle$.

```math
\frac{1}{2} \int \psi_n^\ast(x) V(x) \psi_n(x) \mathrm{d}x = E_n
```
```
 n |        analytical |         numerical 
-- | ----------------- | ----------------- 
 1 |   -1.000000000000 |   -1.000000000000 ✔
 2 |   -0.250000000000 |   -0.250000000000 ✔
 3 |   -0.111111111111 |   -0.111111111111 ✔
 4 |   -0.062500000000 |   -0.062500000000 ✔
 5 |   -0.040000000000 |   -0.040000000000 ✔
 6 |   -0.027777777778 |   -0.027777777778 ✔
 7 |   -0.020408163265 |   -0.020408163265 ✔
 8 |   -0.015625000000 |   -0.015625000000 ✔
 9 |   -0.012345679012 |   -0.012345679012 ✔
10 |   -0.010000000000 |   -0.010000000000 ✔
Test Summary:      | Pass  Total  Time
<ψₙ|V|ψₙ> / 2 = Eₙ |   10     10  0.5s
```

#### Normalization & Orthogonality of $\psi_n(r,\theta,\varphi)$

```math
\int \psi_i^\ast(r,\theta,\varphi) \psi_j(r,\theta,\varphi) r^2 \mathrm{d}r \mathrm{d}\theta \mathrm{d}\varphi = \delta_{ij}
```
```
n₁ | n₂ | l₁ | l₂ | m₁ | m₂ |        analytical |         numerical 
-- | -- | -- | -- | -- | -- | ----------------- | ----------------- 
 1 |  1 |  0 |  0 |  0 |  0 |    1.000000000000 |    1.000000000252 ✔
 1 |  2 |  0 |  0 |  0 |  0 |    0.000000000000 |   -0.000000011223 ✔
 1 |  2 |  0 |  1 |  0 | -1 |    0.000000000000 |   -0.000000000000 ✔
 1 |  2 |  0 |  1 |  0 |  0 |    0.000000000000 |   -0.000000000000 ✔
 1 |  2 |  0 |  1 |  0 |  1 |    0.000000000000 |    0.000000000000 ✔
 1 |  3 |  0 |  0 |  0 |  0 |    0.000000000000 |   -0.000000045661 ✔
 1 |  3 |  0 |  1 |  0 | -1 |    0.000000000000 |    0.000000000000 ✔
 1 |  3 |  0 |  1 |  0 |  0 |    0.000000000000 |   -0.000000000000 ✔
 1 |  3 |  0 |  1 |  0 |  1 |    0.000000000000 |   -0.000000000000 ✔
 1 |  3 |  0 |  2 |  0 | -2 |    0.000000000000 |    0.000000000000 ✔
 1 |  3 |  0 |  2 |  0 | -1 |    0.000000000000 |   -0.000000000000 ✔
 1 |  3 |  0 |  2 |  0 |  0 |    0.000000000000 |   -0.000000000000 ✔
 1 |  3 |  0 |  2 |  0 |  1 |    0.000000000000 |    0.000000000000 ✔
 1 |  3 |  0 |  2 |  0 |  2 |    0.000000000000 |    0.000000000000 ✔
 2 |  1 |  0 |  0 |  0 |  0 |    0.000000000000 |   -0.000000011223 ✔
 2 |  1 |  1 |  0 | -1 |  0 |    0.000000000000 |   -0.000000000000 ✔
 2 |  1 |  1 |  0 |  0 |  0 |    0.000000000000 |   -0.000000000000 ✔
 2 |  1 |  1 |  0 |  1 |  0 |    0.000000000000 |    0.000000000000 ✔
 2 |  2 |  0 |  0 |  0 |  0 |    1.000000000000 |    1.000006970517 ✔
 2 |  2 |  0 |  1 |  0 | -1 |    0.000000000000 |    0.000000000000 ✔
 2 |  2 |  0 |  1 |  0 |  0 |    0.000000000000 |    0.000000000000 ✔
 2 |  2 |  0 |  1 |  0 |  1 |    0.000000000000 |   -0.000000000000 ✔
 2 |  2 |  1 |  0 | -1 |  0 |    0.000000000000 |    0.000000000000 ✔
 2 |  2 |  1 |  0 |  0 |  0 |    0.000000000000 |    0.000000000000 ✔
 2 |  2 |  1 |  0 |  1 |  0 |    0.000000000000 |   -0.000000000000 ✔
 2 |  2 |  1 |  1 | -1 | -1 |    1.000000000000 |    1.000002301351 ✔
 2 |  2 |  1 |  1 | -1 |  0 |    0.000000000000 |   -0.000000000000 ✔
 2 |  2 |  1 |  1 | -1 |  1 |    0.000000000000 |   -0.000000000000 ✔
 2 |  2 |  1 |  1 |  0 | -1 |    0.000000000000 |   -0.000000000000 ✔
 2 |  2 |  1 |  1 |  0 |  0 |    1.000000000000 |    1.000002301351 ✔
 2 |  2 |  1 |  1 |  0 |  1 |    0.000000000000 |    0.000000000000 ✔
 2 |  2 |  1 |  1 |  1 | -1 |    0.000000000000 |   -0.000000000000 ✔
 2 |  2 |  1 |  1 |  1 |  0 |    0.000000000000 |    0.000000000000 ✔
 2 |  2 |  1 |  1 |  1 |  1 |    1.000000000000 |    1.000002301351 ✔
 2 |  3 |  0 |  0 |  0 |  0 |    0.000000000000 |    0.000088519421 ✔
 2 |  3 |  0 |  1 |  0 | -1 |    0.000000000000 |   -0.000000000000 ✔
 2 |  3 |  0 |  1 |  0 |  0 |    0.000000000000 |   -0.000000000000 ✔
 2 |  3 |  0 |  1 |  0 |  1 |    0.000000000000 |    0.000000000000 ✔
 2 |  3 |  0 |  2 |  0 | -2 |    0.000000000000 |   -0.000000000000 ✔
 2 |  3 |  0 |  2 |  0 | -1 |    0.000000000000 |    0.000000000000 ✔
 2 |  3 |  0 |  2 |  0 |  0 |    0.000000000000 |   -0.000000000000 ✔
 2 |  3 |  0 |  2 |  0 |  1 |    0.000000000000 |   -0.000000000000 ✔
 2 |  3 |  0 |  2 |  0 |  2 |    0.000000000000 |   -0.000000000000 ✔
 2 |  3 |  1 |  0 | -1 |  0 |    0.000000000000 |    0.000000000000 ✔
 2 |  3 |  1 |  0 |  0 |  0 |    0.000000000000 |   -0.000000000000 ✔
 2 |  3 |  1 |  0 |  1 |  0 |    0.000000000000 |    0.000000000000 ✔
 2 |  3 |  1 |  1 | -1 | -1 |    0.000000000000 |    0.000038730338 ✔
 2 |  3 |  1 |  1 | -1 |  0 |    0.000000000000 |    0.000000000000 ✔
 2 |  3 |  1 |  1 | -1 |  1 |    0.000000000000 |   -0.000000000000 ✔
 2 |  3 |  1 |  1 |  0 | -1 |    0.000000000000 |   -0.000000000000 ✔
 2 |  3 |  1 |  1 |  0 |  0 |    0.000000000000 |    0.000038730338 ✔
 2 |  3 |  1 |  1 |  0 |  1 |    0.000000000000 |    0.000000000000 ✔
 2 |  3 |  1 |  1 |  1 | -1 |    0.000000000000 |   -0.000000000000 ✔
 2 |  3 |  1 |  1 |  1 |  0 |    0.000000000000 |   -0.000000000000 ✔
 2 |  3 |  1 |  1 |  1 |  1 |    0.000000000000 |    0.000038730338 ✔
 2 |  3 |  1 |  2 | -1 | -2 |    0.000000000000 |   -0.000000000000 ✔
 2 |  3 |  1 |  2 | -1 | -1 |    0.000000000000 |   -0.000000000000 ✔
 2 |  3 |  1 |  2 | -1 |  0 |    0.000000000000 |    0.000000000000 ✔
 2 |  3 |  1 |  2 | -1 |  1 |    0.000000000000 |    0.000000000000 ✔
 2 |  3 |  1 |  2 | -1 |  2 |    0.000000000000 |   -0.000000000272 ✔
 2 |  3 |  1 |  2 |  0 | -2 |    0.000000000000 |    0.000000000000 ✔
 2 |  3 |  1 |  2 |  0 | -1 |    0.000000000000 |    0.000000000000 ✔
 2 |  3 |  1 |  2 |  0 |  0 |    0.000000000000 |   -0.000000000000 ✔
 2 |  3 |  1 |  2 |  0 |  1 |    0.000000000000 |   -0.000000000000 ✔
 2 |  3 |  1 |  2 |  0 |  2 |    0.000000000000 |    0.000000000000 ✔
 2 |  3 |  1 |  2 |  1 | -2 |    0.000000000000 |    0.000000000272 ✔
 2 |  3 |  1 |  2 |  1 | -1 |    0.000000000000 |    0.000000000000 ✔
 2 |  3 |  1 |  2 |  1 |  0 |    0.000000000000 |   -0.000000000000 ✔
 2 |  3 |  1 |  2 |  1 |  1 |    0.000000000000 |   -0.000000000000 ✔
 2 |  3 |  1 |  2 |  1 |  2 |    0.000000000000 |    0.000000000000 ✔
 3 |  1 |  0 |  0 |  0 |  0 |    0.000000000000 |   -0.000000045661 ✔
 3 |  1 |  1 |  0 | -1 |  0 |    0.000000000000 |    0.000000000000 ✔
 3 |  1 |  1 |  0 |  0 |  0 |    0.000000000000 |   -0.000000000000 ✔
 3 |  1 |  1 |  0 |  1 |  0 |    0.000000000000 |    0.000000000000 ✔
 3 |  1 |  2 |  0 | -2 |  0 |    0.000000000000 |    0.000000000000 ✔
 3 |  1 |  2 |  0 | -1 |  0 |    0.000000000000 |    0.000000000000 ✔
 3 |  1 |  2 |  0 |  0 |  0 |    0.000000000000 |   -0.000000000000 ✔
 3 |  1 |  2 |  0 |  1 |  0 |    0.000000000000 |   -0.000000000000 ✔
 3 |  1 |  2 |  0 |  2 |  0 |    0.000000000000 |    0.000000000000 ✔
 3 |  2 |  0 |  0 |  0 |  0 |    0.000000000000 |    0.000088519421 ✔
 3 |  2 |  0 |  1 |  0 | -1 |    0.000000000000 |    0.000000000000 ✔
 3 |  2 |  0 |  1 |  0 |  0 |    0.000000000000 |    0.000000000000 ✔
 3 |  2 |  0 |  1 |  0 |  1 |    0.000000000000 |   -0.000000000000 ✔
 3 |  2 |  1 |  0 | -1 |  0 |    0.000000000000 |   -0.000000000000 ✔
 3 |  2 |  1 |  0 |  0 |  0 |    0.000000000000 |   -0.000000000000 ✔
 3 |  2 |  1 |  0 |  1 |  0 |    0.000000000000 |    0.000000000000 ✔
 3 |  2 |  1 |  1 | -1 | -1 |    0.000000000000 |    0.000038730338 ✔
 3 |  2 |  1 |  1 | -1 |  0 |    0.000000000000 |   -0.000000000000 ✔
 3 |  2 |  1 |  1 | -1 |  1 |    0.000000000000 |   -0.000000000000 ✔
 3 |  2 |  1 |  1 |  0 | -1 |    0.000000000000 |    0.000000000000 ✔
 3 |  2 |  1 |  1 |  0 |  0 |    0.000000000000 |    0.000038730338 ✔
 3 |  2 |  1 |  1 |  0 |  1 |    0.000000000000 |   -0.000000000000 ✔
 3 |  2 |  1 |  1 |  1 | -1 |    0.000000000000 |   -0.000000000000 ✔
 3 |  2 |  1 |  1 |  1 |  0 |    0.000000000000 |    0.000000000000 ✔
 3 |  2 |  1 |  1 |  1 |  1 |    0.000000000000 |    0.000038730338 ✔
 3 |  2 |  2 |  0 | -2 |  0 |    0.000000000000 |   -0.000000000000 ✔
 3 |  2 |  2 |  0 | -1 |  0 |    0.000000000000 |    0.000000000000 ✔
 3 |  2 |  2 |  0 |  0 |  0 |    0.000000000000 |   -0.000000000000 ✔
 3 |  2 |  2 |  0 |  1 |  0 |    0.000000000000 |   -0.000000000000 ✔
 3 |  2 |  2 |  0 |  2 |  0 |    0.000000000000 |   -0.000000000000 ✔
 3 |  2 |  2 |  1 | -2 | -1 |    0.000000000000 |   -0.000000000000 ✔
 3 |  2 |  2 |  1 | -2 |  0 |    0.000000000000 |    0.000000000000 ✔
 3 |  2 |  2 |  1 | -2 |  1 |    0.000000000000 |    0.000000000272 ✔
 3 |  2 |  2 |  1 | -1 | -1 |    0.000000000000 |   -0.000000000000 ✔
 3 |  2 |  2 |  1 | -1 |  0 |    0.000000000000 |    0.000000000000 ✔
 3 |  2 |  2 |  1 | -1 |  1 |    0.000000000000 |    0.000000000000 ✔
 3 |  2 |  2 |  1 |  0 | -1 |    0.000000000000 |    0.000000000000 ✔
 3 |  2 |  2 |  1 |  0 |  0 |    0.000000000000 |    0.000000000000 ✔
 3 |  2 |  2 |  1 |  0 |  1 |    0.000000000000 |   -0.000000000000 ✔
 3 |  2 |  2 |  1 |  1 | -1 |    0.000000000000 |    0.000000000000 ✔
 3 |  2 |  2 |  1 |  1 |  0 |    0.000000000000 |   -0.000000000000 ✔
 3 |  2 |  2 |  1 |  1 |  1 |    0.000000000000 |   -0.000000000000 ✔
 3 |  2 |  2 |  1 |  2 | -1 |    0.000000000000 |   -0.000000000272 ✔
 3 |  2 |  2 |  1 |  2 |  0 |    0.000000000000 |    0.000000000000 ✔
 3 |  2 |  2 |  1 |  2 |  1 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  0 |  0 |  0 |  0 |    1.000000000000 |    1.002052594504 ✔
 3 |  3 |  0 |  1 |  0 | -1 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  0 |  1 |  0 |  0 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  0 |  1 |  0 |  1 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  0 |  2 |  0 | -2 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  0 |  2 |  0 | -1 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  0 |  2 |  0 |  0 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  0 |  2 |  0 |  1 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  0 |  2 |  0 |  2 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  1 |  0 | -1 |  0 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  1 |  0 |  0 |  0 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  1 |  0 |  1 |  0 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  1 |  1 | -1 | -1 |    1.000000000000 |    1.001223346388 ✔
 3 |  3 |  1 |  1 | -1 |  0 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  1 |  1 | -1 |  1 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  1 |  1 |  0 | -1 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  1 |  1 |  0 |  0 |    1.000000000000 |    1.001223346388 ✔
 3 |  3 |  1 |  1 |  0 |  1 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  1 |  1 |  1 | -1 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  1 |  1 |  1 |  0 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  1 |  1 |  1 |  1 |    1.000000000000 |    1.001223346388 ✔
 3 |  3 |  1 |  2 | -1 | -2 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  1 |  2 | -1 | -1 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  1 |  2 | -1 |  0 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  1 |  2 | -1 |  1 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  1 |  2 | -1 |  2 |    0.000000000000 |    0.000000000308 ✔
 3 |  3 |  1 |  2 |  0 | -2 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  1 |  2 |  0 | -1 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  1 |  2 |  0 |  0 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  1 |  2 |  0 |  1 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  1 |  2 |  0 |  2 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  1 |  2 |  1 | -2 |    0.000000000000 |   -0.000000000308 ✔
 3 |  3 |  1 |  2 |  1 | -1 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  1 |  2 |  1 |  0 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  1 |  2 |  1 |  1 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  1 |  2 |  1 |  2 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  2 |  0 | -2 |  0 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  2 |  0 | -1 |  0 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  2 |  0 |  0 |  0 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  2 |  0 |  1 |  0 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  2 |  0 |  2 |  0 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  2 |  1 | -2 | -1 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  2 |  1 | -2 |  0 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  2 |  1 | -2 |  1 |    0.000000000000 |   -0.000000000308 ✔
 3 |  3 |  2 |  1 | -1 | -1 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  2 |  1 | -1 |  0 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  2 |  1 | -1 |  1 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  2 |  1 |  0 | -1 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  2 |  1 |  0 |  0 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  2 |  1 |  0 |  1 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  2 |  1 |  1 | -1 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  2 |  1 |  1 |  0 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  2 |  1 |  1 |  1 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  2 |  1 |  2 | -1 |    0.000000000000 |    0.000000000308 ✔
 3 |  3 |  2 |  1 |  2 |  0 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  2 |  1 |  2 |  1 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  2 |  2 | -2 | -2 |    1.000000000000 |    1.000300628566 ✔
 3 |  3 |  2 |  2 | -2 | -1 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  2 |  2 | -2 |  0 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  2 |  2 | -2 |  1 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  2 |  2 | -2 |  2 |    0.000000000000 |    0.000000193779 ✔
 3 |  3 |  2 |  2 | -1 | -2 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  2 |  2 | -1 | -1 |    1.000000000000 |    1.000300628559 ✔
 3 |  3 |  2 |  2 | -1 |  0 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  2 |  2 | -1 |  1 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  2 |  2 | -1 |  2 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  2 |  2 |  0 | -2 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  2 |  2 |  0 | -1 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  2 |  2 |  0 |  0 |    1.000000000000 |    1.000300628572 ✔
 3 |  3 |  2 |  2 |  0 |  1 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  2 |  2 |  0 |  2 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  2 |  2 |  1 | -2 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  2 |  2 |  1 | -1 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  2 |  2 |  1 |  0 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  2 |  2 |  1 |  1 |    1.000000000000 |    1.000300628559 ✔
 3 |  3 |  2 |  2 |  1 |  2 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  2 |  2 |  2 | -2 |    0.000000000000 |    0.000000193779 ✔
 3 |  3 |  2 |  2 |  2 | -1 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  2 |  2 |  2 |  0 |    0.000000000000 |   -0.000000000000 ✔
 3 |  3 |  2 |  2 |  2 |  1 |    0.000000000000 |    0.000000000000 ✔
 3 |  3 |  2 |  2 |  2 |  2 |    1.000000000000 |    1.000300628566 ✔
Test Summary:                       | Pass  Total  Time
<ψₙ₁ₗ₁ₘ₁|ψₙ₂ₗ₂ₘ₂> = δₙ₁ₙ₂δₗ₁ₗ₂δₘ₁ₘ₂ |  196    196  5.5s

```
