export MorsePotential, V, E, nₘₐₓ, ψ, L

# Packages
using SpecialFunctions

@kwdef struct MorsePotential
  # F. M. Fernández, J. Garcia, ChemistrySelect, 6, 9527−9534(2021) https://doi.org/10.1002/slct.202102509
  # CODATA recommended values of the fundamental physical constants: 2018 https://physics.nist.gov/cgi-bin/cuu/Value?mpsme
  rₑ =  1.997193319969992120068298141276
  Dₑ = - 0.5 - (-0.602634619106539878727562156289)
  k = 2*((-1.1026342144949464615+1/2.00) - (-0.602634619106539878727562156289)) / (2.00 - rₑ)^2
  µ = 1/(1/1836.15267343 + 1/1836.15267343)
  ℏ = 1.0
end

function V(model::MorsePotential, r)
  rₑ = model.rₑ
  Dₑ = model.Dₑ
  k = model.k
  a = sqrt(k/(2*Dₑ))
  if r<0
    throw(DomainError(r, "r=$r is out of the domain (0≦r)"))
  end
  return Dₑ*( exp(-2*a*(r-rₑ)) -2*exp(-a*(r-rₑ)) )
end

function E(model::MorsePotential; n=0)
  Dₑ = model.Dₑ
  k = model.k
  µ = model.µ
  ℏ = model.ℏ
  ω = sqrt(k/µ)
  χ = ℏ*ω/(4*Dₑ)
  return - Dₑ + ℏ*ω*(n+1/2) - χ*ℏ*ω*(n+1/2)^2
end

function nₘₐₓ(model::MorsePotential)
  Dₑ = model.Dₑ
  k = model.k
  µ = model.µ
  ω = sqrt(k/µ)
  return Int(floor((2*Dₑ - ω)/ω))
end

function ψ(model::MorsePotential, r; n=0)
  if r<0
    throw(DomainError(r, "r=$r is out of the domain (0≦r)"))
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

function L(model::MorsePotential, x; n=0, α=0)
  if isinteger(α)
    return sum(k -> (-1)^(k) * (Int(gamma(α+n+1)) // Int((gamma(α+1+k)*gamma(n-k+1)))) * x^k * 1 // factorial(k), 0:n)
  else
    return sum(k -> (-1)^(k) * (gamma(α+n+1) / (gamma(α+1+k)*gamma(n-k+1))) * x^k / factorial(k), 0:n)
  end
end

@doc raw"""
`MorsePotential(α=1.0, m=1.0, ℏ=1.0)`

``\r_e`` is the equilibrium bond distance, ``D_e`` is the the well depth , ``k`` is the force constant, ``\mu`` is the reduced mass and ``\hbar`` is the reduced Planck constant (Dirac's constant).

""" MorsePotential

@doc raw"""
`V(model::MorsePotential, r)`

```math
    V(r) = D_\mathrm{e} \left( \mathrm{e}^{-2a(r-r_e)} - 2\mathrm{e}^{-a(r-r_e)} \right)
```
""" V(model::MorsePotential, r)

@doc raw"""
`E(model::MorsePotential; n=0)`

```math
E_n = - D_\mathrm{e} + \hbar \omega \left( n + \frac{1}{2} \right) - \chi \hbar \omega \left( n + \frac{1}{2} \right)^2
```
""" E(model::MorsePotential; n=0)

@doc raw"""
`nₘₐₓ(model::MorsePotential)`

```math
n_{max} = \left\lfloor \frac{2 D_e - \omega}{\omega} \right\rfloor
```
""" nₘₐₓ(model::MorsePotential)

@doc raw"""
`ψ(model::MorsePotential, r; n=0)`

```math
\psi_n(r) = N_n z^{\lambda-n-1/2} \mathrm{e}^{-z/2} L_n^{(2\lambda-2n-1)}(\xi)```
""" ψ(model::MorsePotential, r; n=0)

@doc raw"""
`L(model::MorsePotential, x; n=0, α=0)`

```math
L_n^{(\alpha)}(x) = \sum_{k=0}^n(-1)^k \frac{\Gamma(\alpha+n+1)}{\Gamma(\alpha+k+1)\Gamma(n-k+1)} \frac{x^k}{k !}
```
""" L(model::MorsePotential, x; n=0, α=0)