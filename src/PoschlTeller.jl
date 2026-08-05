module PoschlTellers

import ..AbstractModel
import ..energy, ..potential, ..wavefunction, ..n_max, ..legendre_polynomial

export PoschlTeller, energy, potential, wavefunction, n_max, legendre_polynomial

# packages
using SpecialFunctions

# parameters
struct PoschlTeller <: AbstractModel
  lambda::Integer
  m::Real
  hbar::Real
  x_0::Real
end

function PoschlTeller(;
  lambda::Integer=1,
  m::Real=1.0,
  hbar::Real=1.0,
  x_0::Real=1.0,
)
  return PoschlTeller(lambda, m, hbar, x_0)
end

# potential
function potential(model::PoschlTeller, x)
  λ  = model.lambda
  m  = model.m
  ℏ  = model.hbar
  x₀ = model.x_0
  return -ℏ^2/(m*x₀^2) * λ*(λ+1)/2 / cosh(x/x₀)^2
end

# maximum quantum number
function n_max(model::PoschlTeller)
  λ = model.lambda
  return Int(floor(λ-1)) # if counting n from zero
end

# eigenvalues
function energy(model::PoschlTeller; n::Integer=0, nocheck=false)
  max_n = n_max(model)
  if !(0 ≤ n ≤ n_max(model) || nocheck)
    throw(DomainError("(n,n_max(model)) = ($n,$(n_max(model)))", "This function is defined for 0 ≤ n ≤ n_max(model)."))
  end
  λ  = model.lambda
  m  = model.m
  ℏ  = model.hbar
  x₀ = model.x_0
  μ  = max_n - n + 1
  return -ℏ^2/(m*x₀^2) * μ^2/2
end

# eigenfunctions
function wavefunction(model::PoschlTeller, x; n::Integer=0)
  max_n = n_max(model)
  if !(0 ≤ n ≤ n_max(model))
    throw(DomainError("(n,n_max(model)) = ($n,$(n_max(model)))", "This function is defined for 0 ≤ n ≤ n_max(model)."))
  end
  λ  = model.lambda
  x₀ = model.x_0
  μ  = max_n - n + 1
  return (-1)^μ / sqrt(x₀) * legendre_polynomial(model,tanh(x/x₀),n=Int64(λ),m=μ) * sqrt(μ*gamma(λ-μ+1)/gamma(λ+μ+1))
end

# associated Legendre polynomials
function legendre_polynomial(model::PoschlTeller, x; n=0, m=0) # same definition as in hydrogen atom: additional factor (-1)^m taken out
  return (1//2)^n * (1-x^2)^(m//2) * sum((-1)^j * factorial(2*n-2*j) // (factorial(j) * factorial(n-j) * factorial(n-2*j-m)) * x^(n-2*j-m) for j ∈ 0:Int(floor((n-m)/2)))
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
  \hat{H} =
  - \frac{\hbar^2}{2 m} \frac{\mathrm{d}^2}{\mathrm{d}x ^2}
  - \frac{\hbar^2}{m x_0^2} \frac{\lambda(\lambda+1)}{2} \frac{1}{\mathrm{cosh}^2(x/x_0)}.
```

After introducing the dimensionless variables
```math
  x^\ast \equiv x/x_0,\qquad E^\ast \equiv \frac{\hbar^2}{m x_0^2} E
```
the Schrödinger equation reduces to
```math
  \hat{H}^\ast \psi(x^\ast) = E^\ast \psi(x^\ast),
```
with
```math
  \hat{H}^\ast = - \frac{1}{2} \frac{\mathrm{d}^2}{\mathrm{d}{x^\ast}^2} - \frac{\lambda(\lambda+1)}{2} \frac{1}{\mathrm{cosh}^2(x^\ast)}.
```
Parameters are specified within the following struct:

```
PT = PoschlTeller(lambda=1, m=1.0, hbar=1.0, x_0=1.0)
```

``\lambda`` determines the potential strength. Currently only integer values for ``\lambda`` are supported.

## References

* [Poschl1933](@cite) G. Pöschl, E. Teller, _Zeitschrift für Physik_, **83** (3–4), 143 (1933), (https://doi.org/10.1007%2FBF01331132). More general definitions are given as (2a), (2b) or (11)
* [Flugge1999](@cite) S. Flügge, Practical Quantum Mechanics (Springer Berlin Heidelberg, 1999), (https://doi.org/10.1007/978-3-642-61995-3) p.94 Problem 39. Potential hole of modified Poschl-Teller type
""" PoschlTeller

@doc raw"""
`potential(model::PoschlTeller, x)`

```math
\begin{aligned}
V(x)
&= -\frac{\hbar^2}{m x_0^2} \frac{\lambda(\lambda+1)}{2} \mathrm{sech}^2(x/x_0) \\
&= -\frac{\hbar^2}{m x_0^2} \frac{\lambda(\lambda+1)}{2} \frac{1}{\mathrm{cosh}^2(x/x_0)} .
\end{aligned}
```
""" potential(model::PoschlTeller, x)

@doc raw"""
`n_max(model::PoschlTeller)`

!!! note
    Note that the number of bound states `n_max + 1` is not equal to the maximum quantum number `n_max`, since we count the ground state from `n=0` in this model.

```math
n_\mathrm{max} = \left\lfloor \lambda \right\rfloor - 1.
```
""" n_max(model::PoschlTeller)

@doc raw"""
`energy(model::PoschlTeller; n::Integer=0, nocheck=false)`

```math
E_n = -\frac{\hbar^2}{m x_0^2}\frac{\mu^2}{2},
```
where ``\mu = \mu(n) = n_\mathrm{max}-n+1``, and ``n_\mathrm{max} = \left\lfloor \lambda \right\rfloor - 1``.
""" energy(model::PoschlTeller; n::Integer=0, nocheck=false)

@doc raw"""
`wavefunction(model::PoschlTeller, x; n::Integer=0)`

```math
\psi_n(x) = \frac{(-1)^μ}{\sqrt{x_0}} P_\lambda^{\mu}(\mathrm{tanh}(x/x_0)) \sqrt{\mu\frac{\Gamma(\lambda-\mu+1)}{\Gamma(\lambda+\mu+1)}},
```
where ``\mu = \mu(n) = n_\mathrm{max}-n+1``, and ``n_\mathrm{max} = \left\lfloor \lambda \right\rfloor - 1`` and ``P_\lambda^{\mu}`` are the associated Legendre functions.
""" wavefunction(model::PoschlTeller, x; n::Integer=0)

@doc raw"""
`legendre_polynomial(model::PoschlTeller, x; n=0, m=0)`

!!! note
    The closed-form expression has been verified to be numerically stable up to $n=10$. For larger $n$, numerical instabilities may arise; a recurrence-relation implementation will address this limitation.

Associated Legendre polynomials are the associated Legendre functions for integer indices. Here we use the same notation of the associated Legendre functions as in the model HydrogenAtom.

```math
\begin{aligned}
P_n^m(x)
&= \left( 1-x^2 \right)^{m/2} \frac{\mathrm{d}^m}{\mathrm{d}x^m} P_n(x) \\
&= \left( 1-x^2 \right)^{m/2} \frac{\mathrm{d}^m}{\mathrm{d}x^m} \frac{1}{2^n n!} \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left[ \left( x^2-1 \right)^n \right] \\
&= \frac{1}{2^n} (1-x^2)^{m/2} \sum_{j=0}^{\left\lfloor\frac{n-m}{2}\right\rfloor} (-1)^j \frac{(2n-2j)!}{j! (n-j)! (n-2j-m)!} x^{(n-2j-m)}.
\end{aligned}
```

""" legendre_polynomial(model::PoschlTeller, x; n=0, m=0)

end # module PoschlTellers