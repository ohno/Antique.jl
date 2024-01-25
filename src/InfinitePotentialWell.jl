export InfinitePotentialWell, V, E, ψ

# Parameters
@kwdef struct InfinitePotentialWell
  L = 1.0
  m = 1.0
  ℏ = 1.0
end

# Potential
function V(model::InfinitePotentialWell, x)
  L = model.L
  return 0<x<L ? 0 : Inf
end

# Energy
function E(model::InfinitePotentialWell; n=1)
  L = model.L
  m = model.m
  ℏ = model.ℏ
  return (ℏ^2*n^2*π^2) / (2*m*L^2)
end

# Wave Function
function ψ(model::InfinitePotentialWell, x; n=1)
  L = model.L
  return 0<x<L ? sqrt(2/L) * sin(n*π*x/L) : 0
end