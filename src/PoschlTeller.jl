export PoschlTeller, V, E, ψ

# parameters
@kwdef struct PoschlTeller
  λ::Int = 1 # Currently only integer values for λ are supported.
  m = 1.0
  ℏ = 1.0
  x₀ = 1.0
end

# potential
function V(model::PoschlTeller, x)
  λ = model.λ
  return -λ*(λ+1)/2/cosh(x)^2
end

# maximum quantum number
function nₘₐₓ(model::PoschlTeller)
  λ = model.λ
  return Int(floor(λ-1)) # if counting n from zero
end

# eigenvalues
function E(model::PoschlTeller; n::Int=0, nocheck=false)
  n_max = nₘₐₓ(model)
  if !(0 ≤ n ≤ nₘₐₓ(model) || nocheck)
    throw(DomainError("(n,nₘₐₓ(model)) = ($n,$(nₘₐₓ(model)))", "This function is defined for 0 ≤ n ≤ nₘₐₓ(model)."))
  end
  λ = model.λ
  m = model.m
  ℏ = model.ℏ
  x₀ = model.x₀
  mu = n_max - n + 1
  return -(mu)^2/2 * ℏ^2/(m*x₀^2)
end

# eigenfunctions
function ψ(model::PoschlTeller, x; n::Int=0)
  n_max = nₘₐₓ(model)
  if !(0 ≤ n ≤ nₘₐₓ(model))
    throw(DomainError("(n,nₘₐₓ(model)) = ($n,$(nₘₐₓ(model)))", "This function is defined for 0 ≤ n ≤ nₘₐₓ(model)."))
  end
  λ = model.λ
  mu = n_max - n + 1
  return (-1)^mu * P(model,tanh(x),n=Int64(λ),m=mu) * sqrt(mu*gamma(λ-mu+1)/gamma(λ+mu+1))
end

# associated Legendre polynomials
function P(model::PoschlTeller, x; n=0, m=0) # same definition as in hydrogen atom: additional factor (-1)^m taken out
  return (1//2)^n * (1-x^2)^(m//2) * sum(j -> (-1)^j * factorial(2*n-2*j) // (factorial(j) * factorial(n-j) * factorial(n-2*j-m)) * x^(n-2*j-m), 0:Int(floor((n-m)/2)))
end

# docstrings

@doc raw"""
`PoschlTeller(λ=1.0, m=1.0, ℏ=1.0, x₀=1.0)`

``\lambda`` determines the potential strength. Currently only integer values for ``\lambda`` are supported.
""" PoschlTeller

@doc raw"""
`V(model::PoschlTeller, x)`

```math
\begin{aligned}
V(x)
&= -\frac{\hbar^2}{m x_0^2} \frac{\lambda(\lambda+1)}{2}  \mathrm{sech}(x)^2
&= -\frac{\hbar^2}{m x_0^2} \frac{\lambda(\lambda+1)}{2}  \frac{1}{\mathrm{cosh}(x)^2}.
\end{aligned}
```
""" V(model::PoschlTeller, x)


@doc raw"""
`nₘₐₓ(model::PoschlTeller)`

```math
n_\mathrm{max} = \left\lfloor \lambda \right\rfloor - 1.
```
Note that the number of bound states is `nₘₐₓ + 1`, since we count the ground state from `n=0`.
""" nₘₐₓ(model::PoschlTeller)

@doc raw"""
`E(model::PoschlTeller; n::Int=0, nocheck=false)`

```math
E_n = -\frac{\hbar^2}{m x_0^2}\frac{\mu^2}{2},
```
where ``\mu = \mu(n) = n_\mathrm{max}-n+1``, and ``n_\mathrm{max} = \left\lfloor \lambda \right\rfloor - 1``.
""" E(model::PoschlTeller; n::Int=0, nocheck=false)

@doc raw"""
`ψ(model::PoschlTeller, x; n::Int=0)`

```math
\psi_n(x) = P_\lambda^{\mu}(\mathrm{tanh}(x/x_0)) \sqrt{\mu\frac{\Gamma(\lambda-\mu+1)}{\Gamma(\lambda+\mu+1)}},
```
where ``\mu = \mu(n) = n_\mathrm{max}-n+1``, and ``n_\mathrm{max} = \left\lfloor \lambda \right\rfloor - 1`` and ``P_\lambda^{\mu}`` are the associated Legendre functions.
""" ψ(model::PoschlTeller, x; n::Int=0)

@doc raw"""
`P(model::PoschlTeller, x; n=0, m=0)`

Associated Legendre polynomials are the associated Legendre functions for integer indices. Here we use the same notation of the associated Legendre functions as in the model HydrogenAtom.
  
```math
\begin{aligned}
P_n^m(x)
&= \left( 1-x^2 \right)^{m/2} \frac{\mathrm{d}^m}{\mathrm{d}x^m} P_n(x) \\
&=  \left( 1-x^2 \right)^{m/2} \frac{\mathrm{d}^m}{\mathrm{d}x^m} \frac{1}{2^n n!} \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left[ \left( x^2-1 \right)^n \right] \\
&= \frac{1}{2^n} (1-x^2)^{m/2} \sum_{j=0}^{\left\lfloor\frac{n-m}{2}\right\rfloor} (-1)^j \frac{(2n-2j)!}{j! (n-j)! (n-2j-m)!} x^{(n-2j-m)}.
\end{aligned}
```
  
""" P(model::PoschlTeller, x; n=0, m=0)