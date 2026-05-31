module MorsePotentials

import ..AbstractModel
import ..energy, ..potential, ..wavefunction, ..n_max, ..laguerre

export MorsePotential, energy, potential, wavefunction, n_max, laguerre

# packages
using SpecialFunctions

# parameters
struct MorsePotential <: AbstractModel
  r_e::Float64
  D_e::Float64
  k::Float64
  mu::Float64
  hbar::Float64
end

function MorsePotential(;
  # The simplified parameters for H2+
  # F. M. Fernandez, J. Garcia, ChemistrySelect, 6, 9527-9534(2021) https://doi.org/10.1002/slct.202102509
  # CODATA recommended values of the fundamental physical constants: 2018 https://physics.nist.gov/cgi-bin/cuu/Value?mpsme
  r_e=2.0, rₑ=r_e,
  D_e=0.1, Dₑ=D_e,
  k=0.1,
  mu=918.1, μ=mu,
  hbar=1.0, ℏ=hbar,
)
  return MorsePotential(rₑ, Dₑ, k, μ, ℏ)
end

function Base.getproperty(model::MorsePotential, sym::Symbol)
  sym === :rₑ && return getfield(model, :r_e)
  sym === :Dₑ && return getfield(model, :D_e)
  sym === :μ && return getfield(model, :mu)
  sym === :ℏ && return getfield(model, :hbar)
  return getfield(model, sym)
end

# potential
function potential(model::MorsePotential, r)
  if !(0 ≤ r)
    throw(DomainError("r = $r", "r must be non-negative: 0 ≤ r."))
  end
  r_e = model.r_e
  D_e = model.D_e
  k = model.k
  a = sqrt(k/(2*D_e))
  return D_e*( exp(-2*a*(r-r_e)) -2*exp(-a*(r-r_e)) )
end

# eigenvalues
function energy(model::MorsePotential; n::Int=0, nocheck=false)
  if !(0 ≤ n ≤ n_max(model) || nocheck)
    throw(DomainError("(n,n_max(model)) = ($n,$(n_max(model)))", "This function is defined for 0 ≤ n ≤ n_max(model)."))
  end
  Dₑ = model.D_e
  k = model.k
  µ = model.mu
  ℏ = model.hbar
  ω = sqrt(k/µ)
  χ = ℏ*ω/(4*Dₑ)
  return - Dₑ + ℏ*ω*(n+1/2) - χ*ℏ*ω*(n+1/2)^2
end

# maximum quantum number
function n_max(model::MorsePotential)
  Dₑ = model.D_e
  k = model.k
  µ = model.mu
  ω = sqrt(k/µ)
  return Int(floor((2*Dₑ - ω)/ω))
end

# eigenfunctions
function wavefunction(model::MorsePotential, r; n::Int=0)
  if !(0 ≤ n ≤ n_max(model))
    throw(DomainError("(n,n_max(model)) = ($n,$(n_max(model)))", "This function is defined for 0 ≤ n ≤ n_max(model)."))
  end
  if !(0 ≤ r)
    throw(DomainError("r = $r", "r must be non-negative: 0 ≤ r."))
  end
  rₑ = model.r_e
  Dₑ = model.D_e
  k = model.k
  µ = model.mu
  ℏ = model.hbar
  a = sqrt(k/(2*Dₑ))
  λ = sqrt(2*µ*Dₑ) / (a*ℏ)
  ξ = 2*λ*exp(-a*(r-rₑ))
  s  = 2*λ - 2*n - 1
  N  = sqrt(factorial(n) * s * a / gamma(s+n+1))
  return N * ξ^(s/2) * exp(-ξ/2) * laguerre(model, ξ, n=n, α=s)
end

# generalized Laguerre polynomials
function laguerre(model::MorsePotential, x; n=0, α=0)
  return laguerre(model, n, α, x)
end
function laguerre(model::MorsePotential, n::Int, α::Int, x)
  return sum((-1)^(k) * (Int(gamma(α+n+1)) // Int((gamma(α+1+k)*gamma(n-k+1)))) * x^k * 1 // factorial(k) for k ∈ 0:n)
end
function laguerre(model::MorsePotential, n::Int, α::Real, x)
  return sum((-1)^(k) * (gamma(α+n+1) / (gamma(α+1+k)*gamma(n-k+1))) * x^k / factorial(k) for k ∈ 0:n)
end

# docstrings

@doc raw"""
## Model

This model is described with the time-independent Schrödinger equation
```math
  \hat{H} \psi(r) = E \psi(r),
```
and the Hamiltonian
```math
  \hat{H} = - \frac{\hbar^2}{2\mu} \frac{\mathrm{d}^2}{\mathrm{d}r ^2} + D_\mathrm{e} \left( \mathrm{e}^{-2a(r-r_e)} - 2\mathrm{e}^{-a(r-r_e)} \right),
```
where ``a = \sqrt{\frac{k}{2Dₑ}}`` is defined. Parameters are specified with the following struct:

```
MP = MorsePotential(rₑ=2.0, Dₑ=0.1, k=0.1, μ=918.1, ℏ=1.0)
```

``r_\mathrm{e}`` is the equilibrium bond distance, ``D_\mathrm{e}`` is the the well depth , ``k`` is the force constant, ``\mu`` is the reduced mass and ``\hbar`` is the reduced Planck constant (Dirac's constant).

## References

* [Morse1929](@cite) P. M. Morse, _Phys. Rev._, **34**, 57 (1929), (https://doi.org/10.1103/PhysRev.34.57)
* [Dahl1988](@cite) J. P. Dahl, M. Springborg, _J. Chem. Phys._, **88**, 4535 (1988). (62), (63) (https://doi.org/10.1063/1.453761)
* [Shao2016](@cite) W. K. Shao, Y. He, J. Pan, _J. Nonlinear Sci. Appl._, **9**, 5, 3388 (2016). (1.6) (http://dx.doi.org/10.22436/jnsa.009.05.124)
* [DLMF](@cite) _The Digital Library of Mathematical Functions_ (DLMF) [18.3 Table1](https://dlmf.nist.gov/18.3#T1), [18.5 Table1](https://dlmf.nist.gov/18.5#T1), [18.5.12](https://dlmf.nist.gov/18.5#E12), [18.5.17_5](https://dlmf.nist.gov/18.5#E17_5)
""" MorsePotential

@doc raw"""
`potential(model::MorsePotential, r)`

```math
V(r) = D_\mathrm{e} \left( \mathrm{e}^{-2a(r-r_e)} - 2\mathrm{e}^{-a(r-r_e)} \right),
```
where ``a = \sqrt{\frac{k}{2Dₑ}}`` is defined. The domain is $0\leq r \lt \infty$.
""" potential(model::MorsePotential, r)

@doc raw"""
`energy(model::MorsePotential; n::Int=0, nocheck=false)`

```math
E_n = - D_\mathrm{e} + \hbar \omega \left( n + \frac{1}{2} \right) - \chi \hbar \omega \left( n + \frac{1}{2} \right)^2,
```
where ``\omega = \sqrt{k/µ}`` and ``\chi = \frac{\hbar\omega}{4D_\mathrm{e}}`` are defined.
""" energy(model::MorsePotential; n::Int=0, nocheck=false)

@doc raw"""
`n_max(model::MorsePotential)`

!!! note
    Note that the number of bound states is equal to the maximum quantum number `n_max`, since we count the ground state from `n=1` in this model.

```math
n_\mathrm{max} = \left\lfloor \frac{2 D_e - \omega}{\omega} \right\rfloor,
```
where ``\omega = \sqrt{k/µ}`` is defined.
""" n_max(model::MorsePotential)

@doc raw"""
`wavefunction(model::MorsePotential, r; n::Int=0)`

```math
\psi_n(r) = N_n z^{\lambda-n-1/2} \mathrm{e}^{-z/2} L_n^{(2\lambda-2n-1)}(\xi),
```

``N_n = \sqrt{\frac{n!(2\lambda-2n-1)a}{\Gamma(2\lambda-n)}}``,
``\lambda = \frac{\sqrt{2\mu D_\mathrm{e}}}{a\hbar}``, ``a = \sqrt{\frac{k}{2Dₑ}}``, ``L_n^{(\alpha)}(x) = \frac{x^{-\alpha} \mathrm{e}^x}{n !} \frac{\mathrm{d}^n}{\mathrm{d} x^n}\left(\mathrm{e}^{-x} x^{n+\alpha}\right)``, ``\xi := 2\lambda\mathrm{e}^{-a(r-r_e)}`` are defined. The domain is $0\leq r \lt \infty$.
""" wavefunction(model::MorsePotential, r; n::Int=0)

@doc raw"""
`laguerre(model::MorsePotential, x; n=0, α=0)`

!!! note
    The generalized Laguerre polynomials $L_n^{(\alpha)}(x)$, not the associated Laguerre polynomials $L_n^{k}(x)$, are used in this model.

Rodrigues' formula & closed-form:
```math
\begin{aligned}
  L_n^{(\alpha)}(x)
  &= \frac{x^{-\alpha}e^x}{n!} \frac{d^n}{dx^n}\left(x^{n+\alpha}e^{-x}\right) \\
  &= \sum_{k=0}^n(-1)^k \left(\begin{array}{l} n+\alpha \\ n-k \end{array}\right) \frac{x^k}{k !} \\
  &= \sum_{k=0}^n(-1)^k \frac{\Gamma(\alpha+n+1)}{\Gamma(\alpha+k+1)\Gamma(n-k+1)} \frac{x^k}{k !}.
\end{aligned}
```
Examples:
```math
\begin{aligned}
  L_0^{(0)}(x) &= 1, \\
  L_1^{(0)}(x) &= 1 - x, \\
  L_1^{(1)}(x) &= 2 - x, \\
  L_2^{(0)}(x) &= 1 - 2 x + 1/2 x^{2}, \\
  L_2^{(1)}(x) &= 3 - 3 x + 1/2 x^{2}, \\
  L_2^{(2)}(x) &= 6 - 4 x + 1/2 x^{2}, \\
  L_3^{(0)}(x) &= 1 - 3 x + 3/2 x^{2} - 1/6 x^{3}, \\
  L_3^{(1)}(x) &= 4 - 6 x + 2 x^{2} - 1/6 x^{3}, \\
  L_3^{(2)}(x) &= 10 - 10 x + 5/2 x^{2} - 1/6 x^{3}, \\
  L_3^{(3)}(x) &= 20 - 15 x + 3 x^{2} - 1/6 x^{3}, \\
  L_4^{(0)}(x) &= 1 - 4 x + 3 x^{2} - 2/3 x^{3} + 1/24 x^{4}, \\
  L_4^{(1)}(x) &= 5 - 10 x + 5 x^{2} - 5/6 x^{3} + 1/24 x^{4}, \\
  L_4^{(2)}(x) &= 15 - 20 x + 15/2 x^{2} - 1 x^{3} + 1/24 x^{4}, \\
  L_4^{(3)}(x) &= 35 - 35 x + 21/2 x^{2} - 7/6 x^{3} + 1/24 x^{4}, \\
  L_4^{(4)}(x) &= 70 - 56 x + 14 x^{2} - 4/3 x^{3} + 1/24 x^{4}, \\
  \vdots
\end{aligned}
```
""" laguerre(model::MorsePotential, x; n=0, α=0)

end # module MorsePotentials