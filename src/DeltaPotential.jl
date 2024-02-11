export DeltaPotential, V, E, ψ

@kwdef struct DeltaPotential
  α = 1.0
  m = 1.0
  ℏ = 1.0
end

function V(model::DeltaPotential, x)
  return x==0 ? -Inf : 0
end

function E(model::DeltaPotential)
  α = model.α
  m = model.m
  ℏ = model.ℏ
  return -(m*α^2)/(2*ℏ^2)
end

function ψ(model::DeltaPotential, x)
  α = model.α
  m = model.m
  ℏ = model.ℏ
  return sqrt(m*α)/ℏ * exp.(-m*α*abs.(x)/ℏ^2)
end

@doc raw"""
`DeltaPotential(α=1.0, m=1.0, ℏ=1.0)`

``\alpha`` is the potential strength, ``m`` is the mass of particle and ``\hbar`` is the reduced Planck constant (Dirac's constant).
""" DeltaPotential

@doc raw"""
`V(model::DeltaPotential, x)`

```math
V(x) = -\alpha \delta(x).
```
""" V(model::DeltaPotential, x) 

@doc raw"""
`E(model::DeltaPotential)`

```math
E = - \frac{m\alpha^2}{2\hbar^2}
```
""" E(model::DeltaPotential)

@doc raw"""
`ψ(model::DeltaPotential, x)`

```math
\psi(x) = \frac{\sqrt{m\alpha}}{\hbar} \mathrm{e}^{-m\alpha |x|/\hbar^2}
```
""" ψ(model::DeltaPotential, x)