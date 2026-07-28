module CoulombTwoBodies

import ..AbstractModel
import ..energy, ..potential, ..wavefunction, ..radial_function, ..laguerre_polynomial, ..spherical_harmonic, ..legendre_polynomial

export CoulombTwoBody, energy, potential, wavefunction, radial_function, laguerre_polynomial, spherical_harmonic, legendre_polynomial

# parameters
struct CoulombTwoBody <: AbstractModel
  z_1::Integer
  z_2::Integer
  m_1::Real
  m_2::Real
  m_e::Real
  a_0::Real
  E_h::Real
  hbar::Real
end

function CoulombTwoBody(;
  z_1::Integer=-1, z₁=z_1,
  z_2::Integer=1, z₂=z_2,
  m_1=1.0, m₁=m_1,
  m_2=1.0, m₂=m_2,
  m_e=1.0, mₑ=m_e,
  a_0=1.0, a₀=a_0,
  E_h=1.0, Eₕ=E_h,
  hbar=1.0, ℏ=hbar,
)
  return CoulombTwoBody(z₁, z₂, m₁, m₂, mₑ, a₀, Eₕ, ℏ)
end

function Base.getproperty(model::CoulombTwoBody, sym::Symbol)
  sym === :z₁ && return getfield(model, :z_1)
  sym === :z₂ && return getfield(model, :z_2)
  sym === :m₁ && return getfield(model, :m_1)
  sym === :m₂ && return getfield(model, :m_2)
  sym === :mₑ && return getfield(model, :m_e)
  sym === :a₀ && return getfield(model, :a_0)
  sym === :Eₕ && return getfield(model, :E_h)
  sym === :ℏ && return getfield(model, :hbar)
  return getfield(model, sym)
end

# potential
function potential(model::CoulombTwoBody, r)
  z_1 = model.z_1
  z_2 = model.z_2
  a_0 = model.a_0
  E_h = model.E_h
  if !(0 ≤ r)
    throw(DomainError("r = $r", "r must be non-negative: 0 ≤ r."))
  end
  return z_1*z_2/abs(r/a_0) * E_h
end

# eigenvalues
function energy(model::CoulombTwoBody; n::Integer=1)
  z_1 = model.z_1
  z_2 = model.z_2
  m_1 = model.m_1
  m_2 = model.m_2
  m_e = model.m_e
  μ = (1/m_1 + 1/m_2)^(-1)
  E_h = model.E_h
  if !(1 ≤ n)
    throw(DomainError("n = $n", "n must be 1 or more: 1 ≤ n."))
  end
  if !(z_1*z_2 < 0)
    throw(DomainError("(z_1,z_2) = ($z_1,$z_2)", "This function is defined for z_1*z_2 < 0."))
  end
  if !(0 < m_1 && 0 < m_2)
    throw(DomainError("(m_1,m_2) = ($m_1,$m_2)", "This function is defined for 0 < m_1, 0 < m_2."))
  end
  return -(z_1*z_2)^2/(2*n^2) * μ/m_e * E_h
end

# eigenfunctions
function wavefunction(model::CoulombTwoBody, r, θ, φ; n::Integer=1, l::Integer=0, m::Integer=0)
  z_1 = model.z_1
  z_2 = model.z_2
  m_1 = model.m_1
  m_2 = model.m_2
  if !(1 ≤ n && 0 ≤ l < n && -l ≤ m ≤ l)
    throw(DomainError("(n,l,m) = ($n,$l,$m)", "This function is defined for 1 ≤ n, 0 ≤ l < n and -l ≤ m ≤ l."))
  end
  if !(0 ≤ r && 0 ≤ θ < π && 0 ≤ φ < 2π)
    throw(DomainError("(r,θ,φ) = ($r,$θ,$φ)", "This function is defined for 0 ≤ r, 0 ≤ θ < π, 0 ≤ φ < 2π."))
  end
  if !(z_1*z_2 < 0)
    throw(DomainError("(z_1,z_2) = ($z_1,$z_2)", "This function is defined for z_1*z_2 < 0."))
  end
  if !(0 < m_1 && 0 < m_2)
    throw(DomainError("(m_1,m_2) = ($m_1,$m_2)", "This function is defined for 0 < m_1, 0 < m_2."))
  end
  return radial_function(model, r, n=n, l=l) * spherical_harmonic(model, θ, φ, l=l, m=m)
end

# radial function
function radial_function(model::CoulombTwoBody, r; n=1, l=0)
  z_1 = model.z_1
  z_2 = model.z_2
  a_0 = model.a_0
  m_1 = model.m_1
  m_2 = model.m_2
  m_e = model.m_e
  μ  = (1/m_1 + 1/m_2)^(-1)
  aμ = a_0 * m_e / μ
  Z  = -z_1*z_2
  ρ = 2*Z*r/(n*aμ)
  N = -sqrt( factorial(n-l-1)/(2*n*factorial(n+l)) * (2*Z/(n*aμ))^3 )
  return N*ρ^l * exp(-ρ/2) * laguerre_polynomial(model, ρ, n=n+l, k=2*l+1)
end

# associated Laguerre polynomials
function laguerre_polynomial(model::CoulombTwoBody, x; n=0, k=0)
  return sum((-1)^(m+k) * factorial(n) // (factorial(m) * factorial(m+k) * factorial(n-m-k)) * x^m for m ∈ 0:n-k)
end

# spherical harmonics
function spherical_harmonic(model::CoulombTwoBody, θ, φ; l=0, m=0)
  N = (-1)^((abs(m)+m)/2) * sqrt( (2*l+1)*factorial(l-Int(abs(m))) / (2*factorial(l+Int(abs(m)))) )
  return N * legendre_polynomial(model,cos(θ), n=l, m=Int(abs(m))) * exp(im*m*φ) / sqrt(2*π)
end

# associated Legendre polynomials
function legendre_polynomial(model::CoulombTwoBody, x; n=0, m=0)
  return (1//2)^n * (1-x^2)^(m//2) * sum((-1)^j * factorial(2*n-2*j) // (factorial(j) * factorial(n-j) * factorial(n-2*j-m)) * x^(n-2*j-m) for j ∈ 0:Int(floor((n-m)/2)))
end

# docstrings

@doc raw"""
## Model

This model is described with the time-independent Schrödinger equation
```math
  \hat{H} \psi(\pmb{r}) = E \psi(\pmb{r}),
```
and the Hamiltonian
```math
  \hat{H} = - \frac{\hbar^2}{2\mu} \nabla^2 + \frac{z_1 z_2}{r/a_0} E_\mathrm{h},
```
where $\mu=\left(\frac{1}{m_1}+\frac{1}{m_2}\right)^{-1}$ is the reduced mass of particle 1 and particle 2. The potential includes only Coulomb interaction and it does not include fine or hyperfine interactions in this model. Parameters are specified with the following struct:

```
CTB = CoulombTwoBody(z_1=-1, z_2=1, m_1=1.0, m_2=1.0, m_e=1.0, a_0=1.0, E_h=1.0, hbar=1.0)
```

``z₁`` is the charge number of particle 1, 
``z₂`` is the charge number of particle 2, 
``m₁`` is the mass of particle 1, 
``m₂`` is the mass of particle 2,
``m_\mathrm{e}`` is the electron mass (use the same unit as ``m₁`` and ``m₂``. For example of hydrogen atom, use ``m_\mathrm{e}=9.1093837139\times10^{-31}\mathrm{kg}``, ``m_1=9.1093837139\times10^{-31}\mathrm{kg}`` and ``m_2=1.67262192595\times10^{-27}\mathrm{kg}`` in the IS unit system, use ``~m_\mathrm{e}=1.0~m_\mathrm{e}``, ``m_1=1.0~m_\mathrm{e}`` and ``m_2=1836.152673426~m_\mathrm{e}`` in the atomic unit.),
``a_0`` is the Bohr radius,
``E_\mathrm{h}`` is the Hartree energy and
``\hbar`` is the reduced Planck constant (Dirac's constant).

## References

* [DLMF](@cite) _The Digital Library of Mathematical Functions_ (DLMF), [18.3 Table1](https://dlmf.nist.gov/18.3#T1), [18.5 Table1](https://dlmf.nist.gov/18.5#T1), [18.5.16](https://dlmf.nist.gov/18.5#E16), [18.5.17](https://dlmf.nist.gov/18.5#E17)
* [cpprefjp](@cite) C++ Japanese Reference (C++日本語リファレンス) _cpprefjp_, [assoc_legendre](https://cpprefjp.github.io/reference/cmath/assoc_legendre.html), [assoc_laguerre](https://cpprefjp.github.io/reference/cmath/assoc_laguerre.html)
* [Messiah1961](@cite) A. Messiah, _Quantum Mechanics_ **VOLUME Ⅰ** (North-Holland Publishing Company, 1961), p.412 I. THE HYDROGEN ATOM
* [Griffiths2018](@cite) D. J. Griffiths, D. F. Schroeter, _Introduction to Quantum Mechanics_ **Third Edition** (Cambridge University Press, 2018) (https://doi.org/10.1017/9781316995433) p.143 4.2 THE HYDROGEN ATOM, p.200 Problem 5.1, p.200 Problem 5.2
* [Greiner2001](@cite) W. Greiner, _Quantum Mechanics: An Introduction_ **Forth Edition** (Springer, 2001) (https://doi.org/10.1007/978-3-642-56826-8) p.217 The Hydrogen Atom, p.236 9.5 Spectrum of a Diatomic Molecule
""" CoulombTwoBody

@doc raw"""
`potential(model::CoulombTwoBody, r)`

```math
\begin{aligned}
  V(r)
  &= - \frac{z_1 z_2 e^2}{4\pi\varepsilon_0 r} 
  &= - \frac{e^2}{4\pi\varepsilon_0 a_0} \frac{z_1 z_2}{r/a_0}
  &= - \frac{z_1 z_2}{r/a_0} E_\mathrm{h},
\end{aligned}
```
where ``E_\mathrm{h} = \frac{\hbar^2}{m_\mathrm{e}{a_0}^2} = \frac{e^2}{4\pi\varepsilon_0a_0} = \frac{m_\mathrm{e}e^4}{\left(4\pi\varepsilon_0\right)^2\hbar^2}`` is the Hartree energy, one of atomic unit. The domain is $0\leq r \lt \infty$.
""" potential(model::CoulombTwoBody, r)

@doc raw"""
`energy(model::CoulombTwoBody; n::Integer=1)`

```math
E_n
= -\frac{(z_1 z_2)^2}{2n^2} \frac{\mu}{m_\mathrm{e}} E_\mathrm{h},
```
where $\mu=\left(\frac{1}{m_1}+\frac{1}{m_2}\right)^{-1}$ is the reduced mass of particle 1 and particle 2, ``E_\mathrm{h} = \frac{\hbar^2}{m_\mathrm{e}{a_0}^2} = \frac{e^2}{4\pi\varepsilon_0a_0} = \frac{m_\mathrm{e}e^4}{\left(4\pi\varepsilon_0\right)^2\hbar^2}`` is the Hartree energy, one of atomic unit. About atomic units, see section 3.9.2 of the [IUPAC GreenBook](https://iupac.org/what-we-do/books/greenbook/). In other units, ``E_\mathrm{h} = 27.211~386~245~988(53)~\mathrm{eV}`` from [here](https://physics.nist.gov/cgi-bin/cuu/Value?hrev).
""" energy(model::CoulombTwoBody; n::Integer=1)

@doc raw"""
`wavefunction(model::CoulombTwoBody, r, θ, φ; n::Integer=1, l::Integer=0, m::Integer=0)`

```math
\psi_{nlm}(\pmb{r}) = R_{nl}(r) Y_{lm}(\theta,\varphi)
```
The domain is $0\leq r \lt \infty, 0\leq \theta \lt \pi, 0\leq \varphi \lt 2\pi$.
""" wavefunction(model::CoulombTwoBody, r, θ, φ; n::Integer=1, l::Integer=0, m::Integer=0)

@doc raw"""
`radial_function(model::CoulombTwoBody, r; n=1, l=0)`

```math
R_{nl}(r) = -\sqrt{\frac{(n-l-1)!}{2n(n+l)!} \left(\frac{2Z}{n a_\mu}\right)^3} \left(\frac{2Zr}{n a_\mu}\right)^l \exp \left(-\frac{Zr}{n a_\mu}\right) L_{n+l}^{2l+1} \left(\frac{2Zr}{n a_\mu}\right),
```
where ``\frac{1}{\mu} = \frac{1}{m_1}+\frac{1}{m_2}``, ``a_\mu = a_0 \frac{m_\mathrm{e}}{\mu}``, ``Z = - z_1 z_2``, the Laguerre polynomials are defined as ``L_n(x) = \frac{1}{n!} \mathrm{e}^x \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left( \mathrm{e}^{-x} x^n \right)``, and the associated Laguerre polynomials are defined as ``L_n^{k}(x) = \frac{\mathrm{d}^k}{\mathrm{d}x^k} L_n(x)``. Note that replace ``2n(n+l)!`` with ``2n[(n+l)!]^3`` if the Laguerre polynomials are defined as ``L_n(x) = \mathrm{e}^x \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left( \mathrm{e}^{-x} x^n \right)``. Note that replace ``L_{n+l}^{2l+1}(x)`` with ``- L_{n-l-1}^{2l+1}(x)`` if the associated Laguerre polynomials are defined as ``L_n^{k}(x) = (-1)^k \frac{\mathrm{d}^k}{\mathrm{d}x^k} L_{n+k}(x)``, which we call the generalized Laguerre polynomials. The domain is $0\leq r \lt \infty$.
""" radial_function(model::CoulombTwoBody, r; n=1, l=0)

@doc raw"""
`laguerre_polynomial(model::CoulombTwoBody, x; n=0, k=0)`

!!! note
    The associated Laguerre polynomials $L_n^{k}(x)$, not the generalized Laguerre polynomials $L_n^{(\alpha)}(x)$, are used in this model.

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
""" laguerre_polynomial(model::CoulombTwoBody, x; n=0, k=0)

@doc raw"""
`spherical_harmonic(model::CoulombTwoBody, θ, φ; l=0, m=0)`

```math
Y_{lm}(\theta,\varphi) = (-1)^{\frac{|m|+m}{2}} \sqrt{\frac{2l+1}{4\pi} \frac{(l-|m|)!}{(l+|m|)!}} P_l^{|m|} (\cos\theta) \mathrm{e}^{im\varphi}.
```
The domain is $0\leq \theta \lt \pi, 0\leq \varphi \lt 2\pi$. Note that some variants are connected by 
```math
i^{|m|+m} \sqrt{\frac{(l-|m|)!}{(l+|m|)!}} P_l^{|m|} = (-1)^{\frac{|m|+m}{2}} \sqrt{\frac{(l-|m|)!}{(l+|m|)!}} P_l^{|m|} = (-1)^m \sqrt{\frac{(l-m)!}{(l+m)!}} P_l^{m}.
```

""" spherical_harmonic(model::CoulombTwoBody, θ, φ; l=0, m=0)

@doc raw"""
`legendre_polynomial(model::CoulombTwoBody, x; n=0, m=0)`

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
""" legendre_polynomial(model::CoulombTwoBody, x; n=0, m=0)

end # module CoulombTwoBodies
