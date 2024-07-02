export MorsePotential, V, E, nₘₐₓ, ψ, L

# packages
using SpecialFunctions

# parameters
@kwdef struct MorsePotential
  # F. M. Fernández, J. Garcia, ChemistrySelect, 6, 9527−9534(2021) https://doi.org/10.1002/slct.202102509
  # CODATA recommended values of the fundamental physical constants: 2018 https://physics.nist.gov/cgi-bin/cuu/Value?mpsme
  rₑ =  1.997193319969992120068298141276
  Dₑ = - 0.5 - (-0.602634619106539878727562156289)
  k = 2*((-1.1026342144949464615+1/2.00) - (-0.602634619106539878727562156289)) / (2.00 - rₑ)^2
  µ = 1/(1/1836.15267343 + 1/1836.15267343)
  ℏ = 1.0
end

# potential
function V(model::MorsePotential, r)
  if !(0 ≤ r)
    throw(DomainError("r = $r", "r must be non-negative: 0 ≤ r."))
  end
  rₑ = model.rₑ
  Dₑ = model.Dₑ
  k = model.k
  a = sqrt(k/(2*Dₑ))
  return Dₑ*( exp(-2*a*(r-rₑ)) -2*exp(-a*(r-rₑ)) )
end

# eigenvalues
function E(model::MorsePotential; n::Int=0, nocheck=false)
  if !(0 ≤ n ≤ nₘₐₓ(model) || nocheck)
    throw(DomainError("(n,nₘₐₓ(model)) = ($n,$(nₘₐₓ(model)))", "This function is defined for 0 ≤ n ≤ nₘₐₓ(model)."))
  end
  Dₑ = model.Dₑ
  k = model.k
  µ = model.µ
  ℏ = model.ℏ
  ω = sqrt(k/µ)
  χ = ℏ*ω/(4*Dₑ)
  return - Dₑ + ℏ*ω*(n+1/2) - χ*ℏ*ω*(n+1/2)^2
end

# maximum quantum number
function nₘₐₓ(model::MorsePotential)
  Dₑ = model.Dₑ
  k = model.k
  µ = model.µ
  ω = sqrt(k/µ)
  return Int(floor((2*Dₑ - ω)/ω))
end

# eigenfunctions
function ψ(model::MorsePotential, r; n::Int=0)
  if !(0 ≤ n ≤ nₘₐₓ(model))
    throw(DomainError("(n,nₘₐₓ(model)) = ($n,$(nₘₐₓ(model)))", "This function is defined for 0 ≤ n ≤ nₘₐₓ(model)."))
  end
  if !(0 ≤ r)
    throw(DomainError("r = $r", "r must be non-negative: 0 ≤ r."))
  end
  rₑ = model.rₑ
  Dₑ = model.Dₑ
  k = model.k
  µ = model.µ
  ℏ = model.ℏ
  a = sqrt(k/(2*Dₑ))
  λ = sqrt(2*µ*Dₑ) / (a*ℏ)
  ξ = 2*λ*exp(-a*(r-rₑ))
  s  = 2*λ - 2*n - 1
  N  = sqrt(factorial(n) * s * a / gamma(s+n+1))
  return N * ξ^(s/2) * exp(-ξ/2) * L(model, ξ, n=n, α=s)
end

# generalized Laguerre polynomials
function L(model::MorsePotential, x; n=0, α=0)
  if isinteger(α)
    return sum(k -> (-1)^(k) * (Int(gamma(α+n+1)) // Int((gamma(α+1+k)*gamma(n-k+1)))) * x^k * 1 // factorial(k), 0:n)
  else
    return sum(k -> (-1)^(k) * (gamma(α+n+1) / (gamma(α+1+k)*gamma(n-k+1))) * x^k / factorial(k), 0:n)
  end
end

# docstrings

@doc raw"""
This model is described with the time-independent Schrödinger equation
```math
  \hat{H} \psi(r) = E \psi(r),
```
and the Hamiltonian
```math
  \hat{H} = - \frac{\hbar^2}{2\mu} \frac{\mathrm{d}^2}{\mathrm{d}r ^2} + D_\mathrm{e} \left( \mathrm{e}^{-2a(r-r_e)} - 2\mathrm{e}^{-a(r-r_e)} \right),
```
where ``a = \sqrt{\frac{k}{2Dₑ}}`` is defined. Parameters are specified with the following struct:

`MP = MorsePotential(rₑ=2.0, Dₑ=0.1, k=0.1, µ=918.1, ℏ=1.0)`

``r_\mathrm{e}`` is the equilibrium bond distance, ``D_\mathrm{e}`` is the the well depth , ``k`` is the force constant, ``\mu`` is the reduced mass and ``\hbar`` is the reduced Planck constant (Dirac's constant).

References:
- [P. M. Morse, _Phys. Rev._, **34**, 57 (1929)](https://doi.org/10.1103/PhysRev.34.57)
- [J. P. Dahl, M. Springborg, _J. Chem. Phys._, **88**, 4535 (1988). (62), (63)](https://doi.org/10.1063/1.453761)
- [W. K. Shao, Y. He, J. Pan, _J. Nonlinear Sci. Appl._, **9**, 5, 3388 (2016). (1.6)](http://dx.doi.org/10.22436/jnsa.009.05.124) 
- The Digital Library of Mathematical Functions (DLMF) [18.3 Table1](https://dlmf.nist.gov/18.3#T1), [18.5 Table1](https://dlmf.nist.gov/18.5#T1), [18.5.12](https://dlmf.nist.gov/18.5#E12), [18.5.17_5](https://dlmf.nist.gov/18.5#E17_5)
""" MorsePotential

@doc raw"""
`V(model::MorsePotential, r)`

```math
V(r) = D_\mathrm{e} \left( \mathrm{e}^{-2a(r-r_e)} - 2\mathrm{e}^{-a(r-r_e)} \right),
```
where ``a = \sqrt{\frac{k}{2Dₑ}}`` is defined. The domain is $0\leq r \lt \infty$.
""" V(model::MorsePotential, r)

@doc raw"""
`E(model::MorsePotential; n::Int=0, nocheck=false)`

```math
E_n = - D_\mathrm{e} + \hbar \omega \left( n + \frac{1}{2} \right) - \chi \hbar \omega \left( n + \frac{1}{2} \right)^2,
```
where ``\omega = \sqrt{k/µ}`` and ``\chi = \frac{\hbar\omega}{4D_\mathrm{e}}`` are defined.
""" E(model::MorsePotential; n::Int=0, nocheck=false)

@doc raw"""
`nₘₐₓ(model::MorsePotential)`

```math
n_\mathrm{max} = \left\lfloor \frac{2 D_e - \omega}{\omega} \right\rfloor,
```
where ``\omega = \sqrt{k/µ}`` is defined.
""" nₘₐₓ(model::MorsePotential)

@doc raw"""
`ψ(model::MorsePotential, r; n::Int=0)`

```math
\psi_n(r) = N_n z^{\lambda-n-1/2} \mathrm{e}^{-z/2} L_n^{(2\lambda-2n-1)}(\xi),
```

``N_n = \sqrt{\frac{n!(2\lambda-2n-1)a}{\Gamma(2\lambda-n)}}``,
``\lambda = \frac{\sqrt{2\mu D_\mathrm{e}}}{a\hbar}``, ``a = \sqrt{\frac{k}{2Dₑ}}``, ``L_n^{(\alpha)}(x) = \frac{x^{-\alpha} \mathrm{e}^x}{n !} \frac{\mathrm{d}^n}{\mathrm{d} x^n}\left(\mathrm{e}^{-x} x^{n+\alpha}\right)``, ``\xi := 2\lambda\mathrm{e}^{-a(r-r_e)}`` are defined. The domain is $0\leq r \lt \infty$.
""" ψ(model::MorsePotential, r; n::Int=0)

@doc raw"""
`L(model::MorsePotential, x; n=0, α=0)`

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
""" L(model::MorsePotential, x; n=0, α=0)