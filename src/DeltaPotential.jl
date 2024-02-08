export DeltaPotential, V, E, ψ

# Parameters
@kwdef struct DeltaPotential
  α = 1.0
  m = 1.0
  ℏ = 1.0
end

# Potential
function V(model::DeltaPotential, x)
  return x==0 ? -Inf : 0
end

# Energy
function E(model::DeltaPotential)
  α = model.α
  m = model.m
  ℏ = model.ℏ
  return -(m*α^2)/(2*ℏ^2)
end

# Wave Function
function ψ(model::DeltaPotential, x)
  α = model.α
  m = model.m
  ℏ = model.ℏ
  return sqrt(m*α)/ℏ * exp.(-m*α*abs.(x)/ℏ^2)
end