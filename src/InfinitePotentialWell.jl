export InfinitePotentialWell, V, E, ψ

@kwdef struct InfinitePotentialWell
  L = 1.0
  m = 1.0
  ℏ = 1.0
end

function V(model::InfinitePotentialWell, x)
  L = model.L
  return 0<x<L ? 0 : Inf
end

function E(model::InfinitePotentialWell; n=1)
  L = model.L
  m = model.m
  ℏ = model.ℏ
  return (ℏ^2*n^2*π^2) / (2*m*L^2)
end

function ψ(model::InfinitePotentialWell, x; n=1)
  L = model.L
  return 0<x<L ? sqrt(2/L) * sin(n*π*x/L) : 0
end

@doc raw"""
`InfinitePotentialWell(L=1.0, m=1.0, ℏ=1.0)`

``L`` is the length of the box, ``m`` is the mass of particle and ``\hbar`` is the reduced Planck constant (Dirac's constant).
""" InfinitePotentialWell

@doc raw"""
`V(model::InfinitePotentialWell; x)`

```math
V(x) =
\left\{
  \begin{array}{ll}
  \infty & x \lt 0, L \lt x \\
  0      & 0 \leq x \leq L
  \end{array}
\right.
```
""" V(model::InfinitePotentialWell, x)

@doc raw"""
`E(model::InfinitePotentialWell; n=1)`

```math
E_n = \frac{\hbar^2 n^2 \pi^2}{2 m L^2}
```
""" E(model::InfinitePotentialWell; n=1)

@doc raw"""
`ψ(model::InfinitePotentialWell, x; n=1)`

```math
\psi_n(x) = \sqrt{\frac{2}{L}} \sin \frac{n\pi x}{L}
```
""" ψ(model::InfinitePotentialWell, x; n=1)