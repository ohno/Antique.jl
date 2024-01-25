export HydrogenAtom, V, E, ψ, R, L, Y, P

# Parameters
@kwdef struct HydrogenAtom
  Z = 1
  ℏ = 1.0
  Eₕ = 1.0
  a₀ = 1.0
  mₑ = 1.0
end

# Potential
function V(model::HydrogenAtom, r; Z=Z, a₀=a₀)
  # if r<0
  #   throw(DomainError(r, "r=$r is out of the domain (0≦r)"))
  # end
  Z = model.Z
  a₀ = model.a₀
  return -Z/abs(r/a₀)
end

# Energy
function E(model::HydrogenAtom; n=1)
  Z = model.Z
  Eₕ = model.Eₕ
  return -Z^2/(2*n^2) * Eₕ
end

# Wave Function
function ψ(model::HydrogenAtom, r, θ, φ; n=1, l=0, m=0)
  # if r<0
  #   throw(DomainError(r, "r=$r is out of the domain (0≦r)"))
  # end
  Z = model.Z
  a₀ = model.a₀
  return R(model, r, n=n, l=l) * Y(model, θ, φ, l=l, m=m)
end

# Radial Function
function R(model::HydrogenAtom, r; n=1, l=0)
  # if r<0
  #   throw(DomainError(r, "r=$r is out of the domain (0≦r)"))
  # end
  Z = model.Z
  a₀ = model.a₀
  ρ = 2*Z*abs(r)/(n*a₀)
  N = -sqrt( factorial(n-l-1)/(2*n*factorial(n+l)) * (2*Z/(n*a₀))^3 )
  return N*ρ^l * exp(-ρ/2) * L(model, ρ, n=n+l, k=2*l+1)
end

# Associated Laguerre Polynomials
function L(model::HydrogenAtom, x; n=0, k=0)
  return sum(m -> (-1)^(m+k) * factorial(n) // (factorial(m) * factorial(m+k) * factorial(n-m-k)) * x^m, 0:n-k)
end

# Spherical Harmonics
function Y(model::HydrogenAtom, θ, φ; l=0, m=0)
  N = (im)^(m+abs(m)) * sqrt( (2*l+1)*factorial(l-Int(abs(m))) / (2*factorial(l+Int(abs(m)))) )
  return N * P(model,cos(θ), n=l, m=Int(abs(m))) * exp(im*m*φ) / sqrt(2*π)
end

# Associated Legendre Polynomials
function P(model::HydrogenAtom, x; n=0, m=0)
  return (1//2)^n * (1-x^2)^(m//2) * sum(j -> (-1)^j * factorial(2*n-2*j) // (factorial(j) * factorial(n-j) * factorial(n-2*j-m)) * x^(n-2*j-m), 0:Int(floor((n-m)/2)))
end