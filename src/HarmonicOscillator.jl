export HarmonicOscillator, V, E, ψ, H

# Parameters
@kwdef struct HarmonicOscillator
  k = 1.0
  m = 1.0
  ℏ = 1.0
end

# Potential
function V(model::HarmonicOscillator, x)
  k = model.k
  return 1//2 * k * x^2
end

# Energy
function E(model::HarmonicOscillator; n=0)
  k = model.k
  m = model.m
  ℏ = model.ℏ
  ω = sqrt(k/m)
  return ℏ * ω * (n+1//2)
end

# Wave Function
function ψ(model::HarmonicOscillator, x; n=0)
  k = model.k
  m = model.m
  ℏ = model.ℏ
  ω = sqrt(k/m)
  A = sqrt(1//(factorial(n)*2^n)*(m*ω/(π*ℏ)))^(1//2)
  ξ = sqrt(m*ω/ℏ) * x
  return A * H(model,ξ,n=n) * exp(-ξ^2/2)
end

# Hermite Polynomials
function H(model::HarmonicOscillator, x; n=0)
  return factorial(n) * sum(i -> (-1)^i // (factorial(i)  * factorial(n-2*i)) * (2*x)^(n-2*i), 0:Int(floor(n/2)))
end
