export HydrogenAtom, V, E, ψ

# parameters
@kwdef struct HydrogenAtom
  Z = 1
  mₑ = 1.0
  a₀ = 1.0
  Eₕ = 1.0
  ℏ = 1.0
end

# potential
function V(model::HydrogenAtom, r)
  if !(0 ≤ r)
    throw(DomainError("r = $r", "r must be non-negative: 0 ≤ r."))
  end
  Z = model.Z
  a₀ = model.a₀
  Eₕ = model.Eₕ
  return -Z/abs(r/a₀) * Eₕ
end

# eigenvalues
function E(model::HydrogenAtom; n::Int=1)
  if !(1 ≤ n)
    throw(DomainError("n = $n", "n must be 1 or more: 1 ≤ n."))
  end
  Z = model.Z
  Eₕ = model.Eₕ
  return -Z^2/(2*n^2) * Eₕ
end

# eigenfunctions
function ψ(model::HydrogenAtom, r, θ, φ; n::Int=1, l::Int=0, m::Int=0)
  if !(1 ≤ n && 0 ≤ l < n && -l ≤ m ≤ l)
    throw(DomainError("(n,l,m) = ($n,$l,$m)", "This function is defined for 1 ≤ n, 0 ≤ l < n and -l ≤ m ≤ l."))
  end
  if !(0 ≤ r && 0 ≤ θ < π && 0 ≤ φ < 2π)
    throw(DomainError("(r,θ,φ) = ($r,$θ,$φ)", "This function is defined for 0 ≤ r, 0 ≤ θ < π, 0 ≤ φ < 2π."))
  end
  return R(model, r, n=n, l=l) * Y(model, θ, φ, l=l, m=m)
end

# radial function
function R(model::HydrogenAtom, r; n=1, l=0)
  Z = model.Z
  a₀ = model.a₀
  ρ = 2*Z*r/(n*a₀)
  N = -sqrt( factorial(n-l-1)/(2*n*factorial(n+l)) * (2*Z/(n*a₀))^3 )
  return N*ρ^l * exp(-ρ/2) * L(model, ρ, n=n+l, k=2*l+1)
end

# associated Laguerre polynomials
function L(model::HydrogenAtom, x; n=0, k=0)
  return sum(m -> (-1)^(m+k) * factorial(n) // (factorial(m) * factorial(m+k) * factorial(n-m-k)) * x^m, 0:n-k)
end      

# spherical harmonics
function Y(model::HydrogenAtom, θ, φ; l=0, m=0)
  N = (-1)^((abs(m)+m)/2) * sqrt( (2*l+1)*factorial(l-Int(abs(m))) / (2*factorial(l+Int(abs(m)))) )
  return N * P(model,cos(θ), n=l, m=Int(abs(m))) * exp(im*m*φ) / sqrt(2*π)
end

# associated Legendre polynomials
function P(model::HydrogenAtom, x; n=0, m=0)
  return (1//2)^n * (1-x^2)^(m//2) * sum(j -> (-1)^j * factorial(2*n-2*j) // (factorial(j) * factorial(n-j) * factorial(n-2*j-m)) * x^(n-2*j-m), 0:Int(floor((n-m)/2)))
end

# docstrings

@doc raw"""
This model is described with the time-independent Schrödinger equation
```math
  \hat{H} \psi(\pmb{r}) = E \psi(\pmb{r}),
```
and the Hamiltonian
```math
  \hat{H} = - \frac{\hbar^2}{2\mu} \nabla^2 + \frac{Z}{r/a_0} E_\mathrm{h},
```
where $\mu=\left(\frac{1}{m_\mathrm{e}}+\frac{1}{m_\mathrm{p}}\right)^{-1}$ is the reduced mass of electron $\mathrm{e}$ and proton $\mathrm{p}$. $\mu = m_\mathrm{e}$ holds in the limit $m_\mathrm{p}\rightarrow\infty$. The potential includes only Coulomb interaction and it does not include fine or hyperfine interactions in this model. Parameters are specified with the following struct:

```
HA = HydrogenAtom(Z=1, mₑ=1.0, a₀=1.0, Eₕ=1.0, ℏ=1.0)
```

``Z`` is the atomic number, ``m_\mathrm{e}`` is the electron mass, ``a_0``is the Bohr radius, ``E_\mathrm{h}`` is the Hartree energy and ``\hbar`` is the reduced Planck constant (Dirac's constant).

Main references:
- _The Digital Library of Mathematical Functions_ (DLMF), [18.3 Table1](https://dlmf.nist.gov/18.3#T1), [18.5 Table1](https://dlmf.nist.gov/18.5#T1), [18.5.16](https://dlmf.nist.gov/18.5#E16), [18.5.17](https://dlmf.nist.gov/18.5#E17)
- _cpprefjp_, [assoc_legendre](https://cpprefjp.github.io/reference/cmath/assoc_legendre.html), [assoc_laguerre](https://cpprefjp.github.io/reference/cmath/assoc_laguerre.html)
- A. Messiah, _Quanfum Mechanics_ **VOLUME Ⅰ** (North-Holland Publishing Company, 1961), p.412 (XI.3), p.419 (XI.18) (XI.18a) (XI.18b), p.483 (B.12), p.493 (B.71) (B.72), p.494 (B.81), p495 (B.93)

Supplemental references:
- cpprefjp, [legendre](https://cpprefjp.github.io/reference/cmath/legendre.html), [assoc_legendre](https://cpprefjp.github.io/reference/cmath/assoc_legendre.html), [laguerre](https://cpprefjp.github.io/reference/cmath/laguerre.html), [assoc_laguerre](https://cpprefjp.github.io/reference/cmath/assoc_laguerre.html)
- The Digital Library of Mathematical Functions (DLMF), [18.3 Table1](https://dlmf.nist.gov/18.3#T1), [18.5 Table1](https://dlmf.nist.gov/18.5#T1), [18.5.16](https://dlmf.nist.gov/18.5#E16), [18.5.17](https://dlmf.nist.gov/18.5#E17), [18.5.12](https://dlmf.nist.gov/18.5#E12)
- L. D. Landau, E. M. Lifshitz, Quantum Mechanics (Pergamon Press, 1965), [p.598 (c.1)](https://archive.org/details/ost-physics-landaulifshitz-quantummechanics/page/n611/mode/2up), [p.598 (c.4)](https://archive.org/details/ost-physics-landaulifshitz-quantummechanics/page/n611/mode/2up), [p.603 (d.13)](https://archive.org/details/ost-physics-landaulifshitz-quantummechanics/page/n615/mode/2up), [p.603 (d.13)](https://archive.org/details/ost-physics-landaulifshitz-quantummechanics/page/n615/mode/2up)
- L. I. Schiff, Quantum Mechanics (McGraw-Hill Book Company, 1968), [p.79 (14.12)](https://archive.org/details/ost-physics-schiff-quantummechanics/page/n95/mode/1up), [p.93 (16.19)](https://archive.org/details/ost-physics-schiff-quantummechanics/page/n109/mode/1up)
- A. Messiah, Quanfum Mechanics (Dover Publications, 1999), [p.493 (B.72)](https://archive.org/details/quantummechanics0000mess/page/491/mode/1up), [p.494 Table](https://archive.org/details/quantummechanics0000mess/page/494/mode/1up), [p.493 (B.72)](https://archive.org/details/quantummechanics0000mess/page/491/mode/1up), [p.483 (B.12)](https://archive.org/details/quantummechanics0000mess/page/483/mode/1up), [p.483 (B.12)](https://archive.org/details/quantummechanics0000mess/page/483/mode/1up)
- W. Greiner, Quantum Mechanics: An Introduction Third Edition (Springer, 1994), [p.83 (4)](https://archive.org/details/quantummechanics0001grei_u4x0/page/83/mode/1up), [p.83 (5)](https://archive.org/details/quantummechanics0001grei_u4x0/page/83/mode/1up), [p.149 (21)](https://archive.org/details/quantummechanics0001grei_u4x0/page/149/mode/1up)
- D. J. Griffiths, Introduction to Quantum Mechanics (Prentice Hall, 1995), [p.126 (4.28)](https://archive.org/details/griffiths-introduction-to-quantum-mechanics/page/126/mode/1up), [p.96 Table3.1](https://archive.org/details/griffiths-introduction-to-quantum-mechanics/page/95/mode/1up), [p.126 (4.27)](https://archive.org/details/griffiths-introduction-to-quantum-mechanics/page/126/mode/1up), [p.139 (4.88)](https://archive.org/details/griffiths-introduction-to-quantum-mechanics/page/139/mode/1up), [p.140 Table4.4](https://archive.org/details/griffiths-introduction-to-quantum-mechanics/page/140/mode/1up), [p.139 (4.87)](https://archive.org/details/griffiths-introduction-to-quantum-mechanics/page/139/mode/1up), [p.140 Table4.5](https://archive.org/details/griffiths-introduction-to-quantum-mechanics/page/140/mode/1up)
- D. A. McQuarrie, J. D. Simon, Physical Chemistry: A Molecular Approach (University Science Books, 1997), [p.195 Table6.1](https://archive.org/details/McQuarrieSimonPhysicalChemistrySolutions/McQuarrie_Simon_Physical_Chemistry1997/page/n218/mode/1up), [p.196 (6.26)](https://archive.org/details/McQuarrieSimonPhysicalChemistrySolutions/McQuarrie_Simon_Physical_Chemistry1997/page/n219/mode/1up), [p.196 Table6.2](https://archive.org/details/McQuarrieSimonPhysicalChemistrySolutions/McQuarrie_Simon_Physical_Chemistry1997/page/n220/mode/1up), [p.207 Table6.4](https://archive.org/details/McQuarrieSimonPhysicalChemistrySolutions/McQuarrie_Simon_Physical_Chemistry1997/page/n230/mode/1up)
- P. W. Atkins, J. De Paula, Atkins' Physical Chemistry, 8th edition (W. H. Freeman, 2008), [p.234](https://archive.org/details/atkinsphysicalch00pwat/page/324/mode/2up?q=Laguerre)
- [J. J. Sakurai, J. Napolitano, Modern Quantum Mechanics Third Edition (Cambridge University Press, 2021)](https://doi.org/10.1017/9781108587280), p.245 Problem 3.30.b, 
""" HydrogenAtom

@doc raw"""
`V(model::HydrogenAtom, r)`

```math
\begin{aligned}
  V(r)
  &= - \frac{Ze^2}{4\pi\varepsilon_0 r} 
  &= - \frac{e^2}{4\pi\varepsilon_0 a_0} \frac{Z}{r/a_0}
  &= - \frac{Z}{r/a_0} E_\mathrm{h},
\end{aligned}
```
where ``E_\mathrm{h} = \frac{\hbar^2}{m_\mathrm{e}{a_0}^2} = \frac{e^2}{4\pi\varepsilon_0a_0} = \frac{m_\mathrm{e}e^4}{\left(4\pi\varepsilon_0\right)^2\hbar^2}`` is the Hartree energy, one of atomic unit. The domain is $0\leq r \lt \infty$.
""" V(model::HydrogenAtom, r)

@doc raw"""
`E(model::HydrogenAtom; n::Int=1)`

```math
E_n
= -\frac{m_\mathrm{e} e^4 Z^2}{2n^2(4\pi\varepsilon_0)^2\hbar^2}
= -\frac{Z^2}{2n^2} E_\mathrm{h},
```
where ``E_\mathrm{h} = \frac{\hbar^2}{m_\mathrm{e}{a_0}^2} = \frac{e^2}{4\pi\varepsilon_0a_0} = \frac{m_\mathrm{e}e^4}{\left(4\pi\varepsilon_0\right)^2\hbar^2}`` is the Hartree energy, one of atomic unit. About atomic units, see section 3.9.2 of the [IUPAC GreenBook](https://iupac.org/what-we-do/books/greenbook/). In other units, ``E_\mathrm{h} = 27.211~386~245~988(53)~\mathrm{eV}`` from [here](https://physics.nist.gov/cgi-bin/cuu/Value?hrev).
""" E(model::HydrogenAtom; n::Int=1)

@doc raw"""
`ψ(model::HydrogenAtom, r, θ, φ; n::Int=1, l::Int=0, m::Int=0)`

```math
\psi_{nlm}(\pmb{r}) = R_{nl}(r) Y_{lm}(\theta,\varphi)
```
The domain is $0\leq r \lt \infty, 0\leq \theta \lt \pi, 0\leq \varphi \lt 2\pi$.
""" ψ(model::HydrogenAtom, r, θ, φ; n::Int=1, l::Int=0, m::Int=0)

@doc raw"""
`R(model::HydrogenAtom, r; n=1, l=0)`

```math
R_{nl}(r) = -\sqrt{\frac{(n-l-1)!}{2n(n+l)!} \left(\frac{2Z}{n a_0}\right)^3} \left(\frac{2Zr}{n a_0}\right)^l \exp \left(-\frac{Zr}{n a_0}\right) L_{n+l}^{2l+1} \left(\frac{2Zr}{n a_0}\right),
```
where Laguerre polynomials are defined as ``L_n(x) = \frac{1}{n!} \mathrm{e}^x \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left( \mathrm{e}^{-x} x^n \right)``, and associated Laguerre polynomials are defined as ``L_n^{k}(x) = \frac{\mathrm{d}^k}{\mathrm{d}x^k} L_n(x)``. Note that replace ``2n(n+l)!`` with ``2n[(n+l)!]^3`` if Laguerre polynomials are defined as ``L_n(x) = \mathrm{e}^x \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left( \mathrm{e}^{-x} x^n \right)``. 
The domain is $0\leq r \lt \infty$.
""" R(model::HydrogenAtom, r; n=1, l=0)

@doc raw"""
`L(model::HydrogenAtom, x; n=0, k=0)`

Rodrigues' formula & closed-form:
```math
\begin{aligned}
L_n^{k}(x)
  &= \frac{\mathrm{d}^k}{\mathrm{d}x^k} L_n(x) \\
  &= \frac{\mathrm{d}^k}{\mathrm{d}x^k} \frac{1}{n!} \mathrm{e}^x \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left( \mathrm{e}^{-x} x^n \right) \\
  &= \sum_{m=0}^{n-k} (-1)^{m+k} \frac{n!}{m!(m+k)!(n-m-k)!} x^m \\
  &= (-1)^k L_{n-k}^{(k)}(x),
\end{aligned}
```
where Laguerre polynomials are defined as ``L_n(x)=\frac{1}{n!}\mathrm{e}^x \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left( \mathrm{e}^{-x} x^n \right)``.

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
  \vdots
\end{aligned}
```
""" L(model::HydrogenAtom, x; n=0, k=0)

@doc raw"""
`Y(model::HydrogenAtom, θ, φ; l=0, m=0)`

```math
Y_{lm}(\theta,\varphi) = (-1)^{\frac{|m|+m}{2}} \sqrt{\frac{2l+1}{4\pi} \frac{(l-|m|)!}{(l+|m|)!}} P_l^{|m|} (\cos\theta) \mathrm{e}^{im\varphi}.
```
The domain is $0\leq \theta \lt \pi, 0\leq \varphi \lt 2\pi$. Note that some variants are connected by 
```math
i^{|m|+m} \sqrt{\frac{(l-|m|)!}{(l+|m|)!}} P_l^{|m|} = (-1)^{\frac{|m|+m}{2}} \sqrt{\frac{(l-|m|)!}{(l+|m|)!}} P_l^{|m|} = (-1)^m \sqrt{\frac{(l-m)!}{(l+m)!}} P_l^{m}.
```

""" Y(model::HydrogenAtom, θ, φ; l=0, m=0)

@doc raw"""
`P(model::HydrogenAtom, x; n=0, m=0)`

Rodrigues' formula & closed-form:
```math
\begin{aligned}
  P_n^m(x)
  &= \left( 1-x^2 \right)^{m/2} \frac{\mathrm{d}^m}{\mathrm{d}x^m} P_n(x) \\
  &= \left( 1-x^2 \right)^{m/2} \frac{\mathrm{d}^m}{\mathrm{d}x^m} \frac{1}{2^n n!} \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left[ \left( x^2-1 \right)^n \right] \\
  &= \frac{1}{2^n} (1-x^2)^{m/2} \sum_{j=0}^{\left\lfloor\frac{n-m}{2}\right\rfloor} (-1)^j \frac{(2n-2j)!}{j! (n-j)! (n-2j-m)!} x^{(n-2j-m)}.
\end{aligned},
```
where Legendre polynomials are defined as ``P_n(x) = \frac{1}{2^n n!} \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left[ \left( x^2-1 \right)^n \right]``. Note that ``P_l^{-m} = (-1)^m \frac{(l-m)!}{(l+m)!} P_l^m`` for ``m<0``. (It is not compatible with ``P_k^m(t) = (-1)^m\left( 1-t^2 \right)^{m/2} \frac{\mathrm{d}^m P_k(t)}{\mathrm{d}t^m}`` caused by ``(-1)^m``.) The specific formulae are given below.

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
  & \vdots
\end{aligned}
```
""" P(model::HydrogenAtom, x; n=0, m=0)