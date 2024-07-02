export DeltaPotential, V, E, ψ

# parameters
@kwdef struct DeltaPotential
  α = 1.0
  m = 1.0
  ℏ = 1.0
end

# potential
function V(model::DeltaPotential, x)
  return x==0 ? -Inf : 0
end

# eigenvalues
function E(model::DeltaPotential)
  α = model.α
  m = model.m
  ℏ = model.ℏ
  return -(m*α^2)/(2*ℏ^2)
end

# eigenfunctions
function ψ(model::DeltaPotential, x)
  α = model.α
  m = model.m
  ℏ = model.ℏ
  return sqrt(m*α)/ℏ * exp.(-m*α*abs.(x)/ℏ^2)
end

# docstrings

@doc raw"""
This model is described with the time-independent Schrödinger equation
```math
  \hat{H} \psi(x) = E \psi(x),
```
and the Hamiltonian
```math
  \hat{H} = - \frac{\hbar^2}{2m} \frac{\mathrm{d}^2}{\mathrm{d}x ^2} - \alpha \delta(x).
```
Parameters are specified with the following struct:

`DeltaPotential(α=1.0, m=1.0, ℏ=1.0)`

``\alpha`` is the potential strength, ``m`` is the mass of particle and ``\hbar`` is the reduced Planck constant (Dirac's constant).

References:
- [D. J. Griffiths, D. F. Schroeter, _Introduction to Quantum Mechanics_ **Third Edition** (Cambridge University Press, 2018)](https://doi.org/10.1017/9781316995433) p.63, 2.5.2 The Delta-Function Well
- [UCSD Physics 130, Quantum Physics](https://quantummechanics.ucsd.edu/ph130a/130_notes/node154.html)
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