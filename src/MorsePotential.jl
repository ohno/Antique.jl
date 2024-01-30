export MorsePotential, V, E, nₘₐₓ, ψ, L

# Packages
using SpecialFunctions

# Parameters
@kwdef struct MorsePotential
  # F. M. Fernández, J. Garcia, ChemistrySelect, 6, 9527−9534(2021) https://doi.org/10.1002/slct.202102509
  # CODATA recommended values of the fundamental physical constants: 2018 https://physics.nist.gov/cgi-bin/cuu/Value?mpsme
  rₑ =  1.997193319969992120068298141276
  Dₑ = - 0.5 - (-0.602634619106539878727562156289)
  k = 2*((-1.1026342144949464615+1/2.00) - (-0.602634619106539878727562156289)) / (2.00 - rₑ)^2
  µ = 1/(1/1836.15267343 + 1/1836.15267343)
  ℏ = 1.0
end

# Potential
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

# Energy
function E(model::MorsePotential; n=0)
  Dₑ = model.Dₑ
  k = model.k
  µ = model.µ
  ℏ = model.ℏ
  ω = sqrt(k/µ)
  χ = ℏ*ω/(4*Dₑ)
  return - Dₑ + ℏ*ω*(n+1/2) - χ*ℏ*ω*(n+1/2)^2
end

# Maximum of n
function nₘₐₓ(model::MorsePotential)
  Dₑ = model.Dₑ
  k = model.k
  µ = model.µ
  ω = sqrt(k/µ)
  return Int(floor((2*Dₑ - ω)/ω))
end

# Wave Function
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

# Generalized Laguerre Polynomials
function L(model::MorsePotential, x; n=0, α=0)
  if isinteger(α)
    return sum(k -> (-1)^(k) * (Int(gamma(α+n+1)) // Int((gamma(α+1+k)*gamma(n-k+1)))) * x^k * 1 // factorial(k), 0:n)
  else
    return sum(k -> (-1)^(k) * (gamma(α+n+1) / (gamma(α+1+k)*gamma(n-k+1))) * x^k / factorial(k), 0:n)
  end
end