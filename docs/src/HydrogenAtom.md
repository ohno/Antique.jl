```@meta
CurrentModule = Antique
```

# Hydrogen Atom

The hydrogen atom is the simplest Coulomb 2-body system.

## Definitions

This model is described with the time-independent Schrödinger equation
```math
  \hat{H} \psi(\pmb{r}) = E \psi(\pmb{r}),
```
and the Hamiltonian
```math
  \hat{H} = - \frac{\hbar^2}{2\mu} \nabla^2 + V(r),
```
where $\mu=\left(\frac{1}{m_\mathrm{e}}+\frac{1}{m_\mathrm{p}}\right)^{-1}$ is the reduced mass of electron $\mathrm{e}$ and proton $\mathrm{p}$. $\mu = m_\mathrm{e}$ holds in the limit $m_\mathrm{p}\rightarrow\infty$. The potential includes only Coulomb interaction and it does not include fine or hyperfine interactions in this model. Parameters are specified with the following struct.

#### Parameters
```@docs; canonical=false
Antique.HydrogenAtom
```

#### Potential
```@docs; canonical=false
Antique.V(::HydrogenAtom, ::Any)
```

#### Eigen Values
```@docs; canonical=false
Antique.E(::HydrogenAtom)
```

#### Eigen Functions
```@docs; canonical=false
Antique.ψ(::HydrogenAtom, ::Any, ::Any, ::Any)
```

#### Radial Functions
```@docs; canonical=false
Antique.R(::HydrogenAtom, ::Any)
```

#### Associated Laguerre Polynomials
```@docs; canonical=false
Antique.L(::HydrogenAtom, ::Any)
```

#### Spherical Harmonics
```@docs; canonical=false
Antique.Y(::HydrogenAtom, ::Any, ::Any)
```

#### Associated Legendre Polynomials
```@docs; canonical=false
Antique.P(::HydrogenAtom, ::Any)
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

[Install Antique.jl](@ref Install) for the first use and run `using Antique` before each use. The energy `E()`, wavefunction `ψ()`, potential `V()` and some other functions are suppoted. In this system, the model is generated by `HydrogenAtom` and several parameters `Z`, `Eₕ`, `mₑ`, `a₀` and `ℏ` are set as optional arguments.

```@example HA
using Antique
H = HydrogenAtom(Z=1, Eₕ=1.0, a₀=1.0, mₑ=1.0, ℏ=1.0)
; #hide
```

Parameters:

```@repl HA
H.Z
H.Eₕ
H.mₑ
H.a₀
H.ℏ
```

Eigen values:

```@repl HA
E(H, n=1)
E(H, n=2)
```

Wave length ($n=2\rightarrow1$, the first line of the Lyman series):

```@example HA
Eₕ2nm⁻¹ = 2.1947463136320e-2 # https://physics.nist.gov/cgi-bin/cuu/CCValue?hrminv
println("ΔE = ", E(H,n=2) - E(H,n=1), " Eₕ")
println("λ  = ", ((E(H,n=2)-E(H,n=1))*Eₕ2nm⁻¹)^-1, " nm")
```

Hyperfine Splitting:

```@example HA
# E. Tiesinga, et al., Rev. Mod. Phys. 93, 025010 (2021) https://doi.org/10.1103/RevModPhys.93.025010
e  = 1.602176634e-19    # C      https://physics.nist.gov/cgi-bin/cuu/Value?e
h  = 6.62607015e-34     # J Hz-1 https://physics.nist.gov/cgi-bin/cuu/Value?h
c  = 299792458          # m s-1  https://physics.nist.gov/cgi-bin/cuu/Value?c
a0 = 5.29177210903e-11  # m      https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0
μ0 = 1.25663706212e-6   # N A-2  https://physics.nist.gov/cgi-bin/cuu/Value?mu0
μB = 9.2740100783e-24   # J T-1  https://physics.nist.gov/cgi-bin/cuu/Value?mub
μN = 5.0507837461e-27   # J T-1  https://physics.nist.gov/cgi-bin/cuu/Value?mun
ge = 2.00231930436256   #        https://physics.nist.gov/cgi-bin/cuu/Value?gem
gp = 5.5856946893       #        https://physics.nist.gov/cgi-bin/cuu/Value?gp

# D. J. Griffiths, Am. J. Phys. 50, 698 (1982) https://doi.org/10.1119/1.12733
δ = abs(ψ(H,0,0,0))^2
ΔE = 2 / 3 * μ0 * μN * μB * gp * ge * δ * a0^(-3)
println("1/π    = ", 1/π)
println("<δ(r)> = ", δ, " a₀⁻³")
println("<δ(r)> = ", δ * a0^(-3), " m⁻³")
println("ΔE = ", ΔE, " J")
println("ν = ΔE/h = ", ΔE / h * 1e-6, " MHz")
println("λ = hc/ΔE = ", h*c/ΔE*100, " cm")
```

Potential energy curve:

```@example HA
using CairoMakie

f = Figure()
ax = Axis(f[1,1], xlabel=L"$r~/~a_0$", ylabel=L"$V(r)~/~E_\mathrm{h}$",  limits=(0.0,15.0,-2.0,0.2))
lines!(ax, 0.1:0.01:20, r -> V(H, r))
f
```

Radial functions:

```@example HA
using CairoMakie
using LaTeXStrings

# setting
f = Figure()
ax = Axis(f[1,1], xlabel=L"$r~/~a_0$", ylabel=L"$r^2|R_{nl}(r)|^2~/~a_0^{-1}$", limits=(0,20,0,0.58))

# plot
ws = []
ls = []
for n in 1:3
  for l in 0:n-1
    w = lines!(
        ax,
        0..20,
        r -> r^2 * R(H,r,n=n,l=l)^2,
        linewidth = 2,
        linestyle = [:solid,:dash,:dot,:dashdot,:dashdotdot][l+1],
        color = n,
        colormap = :tab10,
        colorrange = (1,10)
    )
    push!(ws, w)
    push!(ls, latexstring("n=$n, l=$l"))
  end
end

# legend
axislegend(ax, ws, ls, position=:rt)

f
```

Wave functions (electron density in $n=5,l=2,m=1$):

```@example HA
using Antique
H = HydrogenAtom(Z=1, Eₕ=1.0, a₀=1.0, mₑ=1.0, ℏ=1.0)
loop(x) = x<-1 ? loop(x+2) : (1<x ? loop(x-2) : x)
myacos(x) = acos(loop(x))
r(x,y,z)  = sqrt(x^2+y^2+z^2)
θ(x,y,z) = x^2+y^2<1e-9 ? 0 : myacos(z/r(x,y,z)) 
φ(x,y,z) = y^2<1e-9 ? 0 : sign(y)*myacos(x/sqrt(x^2+y^2))
P(x,y,z) = abs(ψ(H,r(x,y,z),θ(x,y,z),φ(x,y,z),n=5,l=2,m=1))^2

using CairoMakie
f = Figure(size=(500,500), backgroundcolor=:transparent)
a = Axis(f[1,1], aspect=1)
hidespines!(a)
hidedecorations!(a)
heatmap!(a, -40:0.1:40, -40:0.1:40, (y,z) -> P(0,y,z), colorrange=(0.0,0.00001))
f
save("assets/fig/HydrogenAtom.png", f) # hide
; # hide
```

![](assets/fig/HydrogenAtom.png)

## Testing

Unit testing and Integration testing were done using computer algebra system ([Symbolics.jl](https://symbolics.juliasymbolics.org/stable/)) and numerical integration ([QuadGK.jl](https://juliamath.github.io/QuadGK.jl/stable/)). The test script is [here](https://github.com/ohno/Antique.jl/blob/main/test/HydrogenAtom.jl).

```@eval
using Markdown
using Antique
Markdown.parse(Antique.load("../../test/result/HydrogenAtom.log"))
```
