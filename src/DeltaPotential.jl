module DeltaPotentials

# for Julia 1.1
import Base:@kwdef

import ..AbstractModel
import ..energy, ..potential, ..wavefunction

export DeltaPotential, energy, potential, wavefunction

# parameters
@kwdef struct DeltaPotential <: AbstractModel
  alpha::Real = 1.0
  m::Real = 1.0
  hbar::Real = 1.0
end

# potential
function potential(model::DeltaPotential, x)
  return x==0 ? -Inf : 0
end

# eigenvalues
function energy(model::DeltaPotential)
  α = model.alpha
  m = model.m
  ℏ = model.hbar
  return -(m*α^2)/(2*ℏ^2)
end

# eigenfunctions
function wavefunction(model::DeltaPotential, x)
  α = model.alpha
  m = model.m
  ℏ = model.hbar
  return sqrt(m*α)/ℏ * exp.(-m*α*abs.(x)/ℏ^2)
end

# docstrings

@doc raw"""
## Model

This model is described with the time-independent Schrödinger equation
```math
  \hat{H} \psi(x) = E \psi(x),
```
and the Hamiltonian
```math
  \hat{H} = - \frac{\hbar^2}{2m} \frac{\mathrm{d}^2}{\mathrm{d}x ^2} - \alpha \delta(x).
```
Parameters are specified with the following struct:

```
DP = DeltaPotential(alpha=1.0, m=1.0, hbar=1.0)
```

``\alpha`` is the potential strength, ``m`` is the mass of particle and ``\hbar`` is the reduced Planck constant (Dirac's constant).

## References

* [Griffiths2018](@cite) D. J. Griffiths, D. F. Schroeter, _Introduction to Quantum Mechanics_ **Third Edition** (Cambridge University Press, 2018), (https://doi.org/10.1017/9781316995433) p.63, 2.5.2 The Delta-Function Well
* [DeltaFunctionPotential](@cite) UCSD Physics 130, Quantum Physics, (https://quantummechanics.ucsd.edu/ph130a/130_notes/node154.html)
""" DeltaPotential

@doc raw"""
`potential(model::DeltaPotential, x)`

```math
V(x) = -\alpha \delta(x).
```
""" potential(model::DeltaPotential, x) 

@doc raw"""
`energy(model::DeltaPotential)`

```math
E = - \frac{m\alpha^2}{2\hbar^2}
```
""" energy(model::DeltaPotential)

@doc raw"""
`wavefunction(model::DeltaPotential, x)`

```math
\psi(x) = \frac{\sqrt{m\alpha}}{\hbar} \mathrm{e}^{-m\alpha |x|/\hbar^2}
```
""" wavefunction(model::DeltaPotential, x)

end # module DeltaPotentials