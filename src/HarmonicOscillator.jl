module HarmonicOscillators

import ..AbstractModel
import ..energy, ..potential, ..wavefunction, ..laguerre_polynomial

export HarmonicOscillator, energy, potential, wavefunction, laguerre_polynomial

# parameters
struct HarmonicOscillator <: AbstractModel
  k::Real
  m::Real
  hbar::Real
end

function HarmonicOscillator(;
  k=1.0,
  m=1.0,
  hbar=1.0, ℏ=hbar,
)
  return HarmonicOscillator(k, m, ℏ)
end

function Base.getproperty(model::HarmonicOscillator, sym::Symbol)
  sym === :ℏ && return getfield(model, :hbar)
  return getfield(model, sym)
end

# potential
function potential(model::HarmonicOscillator, x)
  k = model.k
  return 1//2 * k * x^2
end

# eigenvalues
function energy(model::HarmonicOscillator; n::Integer=0)
  if !(0 ≤ n)
    throw(DomainError("n = $n", "n must be non-negative: 0 ≤ n."))
  end
  k = model.k
  m = model.m
  hbar = model.hbar
  ω = sqrt(k/m)
  return hbar * ω * (n+1//2)
end

# eigenfunctions
function wavefunction(model::HarmonicOscillator, x; n::Integer=0)
  if !(0 ≤ n)
    throw(DomainError("n = $n", "n must be non-negative: 0 ≤ n."))
  end
  k = model.k
  m = model.m
  hbar = model.hbar
  ω = sqrt(k/m)
  A = sqrt(1//(factorial(n)*2^n)*sqrt(m*ω/(π*hbar)))
  ξ = sqrt(m*ω/hbar) * x
  return A * laguerre_polynomial(model,ξ,n=n) * exp(-ξ^2/2)
end

# Hermite polynomials
function laguerre_polynomial(model::HarmonicOscillator, x; n=0)
  return factorial(n) * sum((-1)^i // (factorial(i)  * factorial(n-2*i)) * (2*x)^(n-2*i) for i ∈ 0:Int(floor(n/2)))
end

# docstrings

@doc raw"""
## Model

This model is described with the time-independent Schrödinger equation
```math
  \hat{H} \psi(x) = E \psi(x),
```
and the Hamiltonian
```math
  \hat{H} = - \frac{\hbar^2}{2m} \frac{\mathrm{d}^2}{\mathrm{d}x ^2} + \frac{1}{2} k x^2.
```
Parameters are specified with the following struct:

```
HO = HarmonicOscillator(k=1.0, m=1.0, hbar=1.0)
```

``k`` is the force constant, ``m`` is the mass of the particle and ``\hbar`` is the reduced Planck constant (Dirac's constant).

## References

Main:
* [DLMF](@cite) _The Digital Library of Mathematical Functions_ (DLMF) [18.5.18](https://dlmf.nist.gov/18.5#E18)
* [cpprefjp](@cite) C++ Japanese Reference _cpprefjp_, [laguerre_polynomial](https://cpprefjp.github.io/reference/cmath/laguerre_polynomial.html)
* [Griffiths2018](@cite) D. J. Griffiths, D. F. Schroeter, _Introduction to Quantum Mechanics_ **Third Edition** (Cambridge University Press, 2018) p.48, 2.3.2 Analytic Method

Supplemental:
* [DLMF](@cite) The Digital Library of Mathematical Functions (DLMF)                                                    [18.3 Table1](https://dlmf.nist.gov/18.3#T1), [18.5 Table1](https://dlmf.nist.gov/18.5#T1), [18.5.13](https://dlmf.nist.gov/18.5#E13), [18.5.18](https://dlmf.nist.gov/18.5#E18)
* [Landau1965](@cite) L. D. Landau, E. M. Lifshitz, Quantum Mechanics (Pergamon Press, 1965)                                  [p.595 (a.4), (a.6)](https://archive.org/details/ost-physics-landaulifshitz-quantummechanics/page/n607/mode/2up)
* [Schiff1968](@cite) L. I. Schiff, Quantum Mechanics (McGraw-Hill Book Company, 1968)                                        [p.71 (13.12)](https://archive.org/details/ost-physics-schiff-quantummechanics/page/n87/mode/2up)
* [Messiah1999](@cite) A. Messiah, Quanfum Mechanics (Dover Publications, 1999)                                                [p.491 (B.59)](https://archive.org/details/quantummechanics0000mess/page/491/mode/1up)
* [Greiner1994](@cite) W. Greiner, Quantum Mechanics: An Introduction Third Edition (Springer, 1994)                           [p.152 (7.22)](https://archive.org/details/quantummechanics0001grei_u4x0/page/152/mode/1up)
* [Griffiths1995](@cite) D. J. Griffiths, Introduction to Quantum Mechanics (Prentice Hall, 1995)                                [p.41 Table 2.1](https://archive.org/details/griffiths-introduction-to-quantum-mechanics/page/41/mode/1up), [p.43 (2.70)](https://archive.org/details/griffiths-introduction-to-quantum-mechanics/page/43/mode/1up)
* [McQuarrie1997](@cite) D. A. McQuarrie, J. D. Simon, Physical Chemistry: A Molecular Approach (University Science Books, 1997) [p.170 Table 5.2](https://archive.org/details/McQuarrieSimonPhysicalChemistrySolutions/McQuarrie_Simon_Physical_Chemistry1997/page/n193/mode/1up)
* [Atkins2008](@cite) P. W. Atkins, J. De Paula, Atkins' Physical Chemistry, 8th edition (W. H. Freeman, 2008)                [p.293 Table 9.1](https://archive.org/details/atkinsphysicalch00pwat/page/292/mode/2up)
* [Sakurai2021](@cite) J. J. Sakurai, J. Napolitano, Modern Quantum Mechanics Third Edition (Cambridge University Press, 2021) [p.524 (B.29)](https://doi.org/10.1017/9781108587280)
""" HarmonicOscillator

@doc raw"""
`potential(model::HarmonicOscillator, x)`

```math
V(x)
= \frac{1}{2} k x^2
= \frac{1}{2} m \omega^2 x^2
= \frac{1}{2} \hbar \omega \xi^2,
```
where ``\omega = \sqrt{k/m}`` is the angular frequency and ``\xi = \sqrt{\frac{m\omega}{\hbar}}x``.
""" potential(model::HarmonicOscillator, x)

@doc raw"""
`energy(model::HarmonicOscillator; n::Integer=0)`

```math
E_n = \hbar \omega \left( n + \frac{1}{2} \right),
```
where ``\omega = \sqrt{k/m}`` is the angular frequency.
""" energy(model::HarmonicOscillator; n::Integer=0)

@doc raw"""
`wavefunction(model::HarmonicOscillator, x; n::Integer=0)`

```math
\psi_n(x) = A_n H_n(\xi) \exp{\left( -\frac{\xi^2}{2} \right)},
```
where ``\omega = \sqrt{k/m}``, ``\xi = \sqrt{\frac{m\omega}{\hbar}}x``, ``A_n = \sqrt{\frac{1}{n! 2^n} \sqrt{\frac{m\omega}{\pi\hbar}}}``, ``H_n(x) = (-1)^n \mathrm{e}^{x^2} \frac{\mathrm{d}^n}{\mathrm{d}x^n} \mathrm{e}^{-x^2}`` are defined.
""" wavefunction(model::HarmonicOscillator, x; n::Integer=0)

@doc raw"""
`laguerre_polynomial(model::HarmonicOscillator, x; n=0)`

!!! note
    The closed-form expression has been verified to be numerically stable up to $n=10$. For larger $n$, numerical instabilities may arise; a recurrence-relation implementation will address this limitation.

Rodrigues' formula & closed-form:
```math
\begin{aligned}
  H_{n}(x)
  &:= (-1)^n \mathrm{e}^{x^2} \frac{\mathrm{d}^n}{\mathrm{d}x^n} \mathrm{e}^{-x^2} \\
  &= n! \sum_{m=0}^{\lfloor n/2 \rfloor} \frac{(-1)^m}{m! (n-2m)!}(2 x)^{n-2m}.
\end{aligned}
```
Examples:
```math
\begin{aligned}
  H_{0}(x)  &= 1, \\
  H_{1}(x)  &= 2 x, \\
  H_{2}(x)  &= -2 + 4 x^{2}, \\
  H_{3}(x)  &= -12 x + 8 x^{3}, \\
  H_{4}(x)  &= 12 - 48 x^{2} + 16 x^{4}, \\
  H_{5}(x)  &= 120 x - 160 x^{3} + 32 x^{5}, \\
  H_{6}(x)  &= -120 + 720 x^{2} - 480 x^{4} + 64 x^{6}, \\
  H_{7}(x)  &= -1680 x + 3360 x^{3} - 1344 x^{5} + 128 x^{7}, \\
  H_{8}(x)  &= 1680 - 13440 x^{2} + 13440 x^{4} - 3584 x^{6} + 256 x^{8}, \\
  H_{9}(x)  &= 30240 x - 80640 x^{3} + 48384 x^{5} - 9216 x^{7} + 512 x^{9}, \\
  &\vdots
\end{aligned}
```
""" laguerre_polynomial(model::HarmonicOscillator, x; n=0)

end # module HarmonicOscillators
