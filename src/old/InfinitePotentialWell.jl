module OldInfinitePotentialWell

  # Default
  L = 1.0 # change here!
  m = 1.0 # change here!
  ℏ = 1.0 # change here!

  # Potential
  V(x; L=L) = 0<x<L ? 0 : Inf

  # Wave Function
  ψ(x; n=1, L=L) = 0<x<L ? sqrt(2/L) * sin(n*π*x/L) : 0

  # Energy
  E(; n=1, L=L, m=m, ℏ=ℏ) = (ℏ^2*n^2*π^2) / (2*m*L^2)

end