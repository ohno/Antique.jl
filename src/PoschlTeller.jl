export PoschlTeller, V, E, ψ

@kwdef struct PoschlTeller
  λ = 1.0
  m = 1.0
  ℏ = 1.0
  x₀ = 1.0
end

function inputchk(model::PoschlTeller,λ,n,n_max)
  if (λ ≈ round(λ)) == false
    @show(λ,round(λ),λ ≈ round(λ))
    error("Error: Currently only integer values for λ are supported.")
  end
  
  if (n ≈ round(n)) == false
    @show(n,round(n),n ≈ round(n))
    error("Error: n is the index of the n-th excited state. It must be an integer.")
  end
  
  if n < 0 || n > n_max
    @show(n,n_max)
    error("Error: n must be non-negative and smaller than or equal n_max: 0 <= n <= n_max")
  end
end


function V(model::PoschlTeller, x)
  λ = model.λ
  return -λ*(λ+1)/2/cosh(x)^2
end

function nₘₐₓ(model::PoschlTeller)
  λ = model.λ
  return Int(floor(λ-1)) # if counting n from zero
end

function E(model::PoschlTeller; n=0)
  λ = model.λ
  m = model.m
  ℏ = model.ℏ
  x₀ = model.x₀
  
  n_max = nₘₐₓ(model)
  mu = n_max-n+1
  inputchk(model,λ,n,n_max)
  return -(mu)^2/2 * ℏ^2/(m*x₀^2)
end

function ψ(model::PoschlTeller, x; n=0)
  λ = model.λ
  n_max = nₘₐₓ(model)
  mu = n_max-n+1
  
  inputchk(model,λ,n,n_max)
  return (-1)^mu * P(model,tanh(x),n=Int64(λ),m=mu)*sqrt(mu*gamma(λ-mu+1)/gamma(λ+mu+1))
end

function P(model::PoschlTeller, x; n=0, m=0) # same definition as in hydrogen atom: additional factor (-1)^m taken out
  return (1//2)^n * (1-x^2)^(m//2) * sum(j -> (-1)^j * factorial(2*n-2*j) // (factorial(j) * factorial(n-j) * factorial(n-2*j-m)) * x^(n-2*j-m), 0:Int(floor((n-m)/2)))
end

# docstrings:
@doc raw"""
`PoschlTeller(λ=1.0, m=1.0, ℏ=1.0, x₀=1.0)`

``\lambda`` determines the potential strength.
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
""" nₘₐₓ(model::PoschlTeller)

@doc raw"""
`E(model::PoschlTeller; n=0)`

```math
E_n = -\frac{\hbar^2}{m x_0^2}\frac{\mu^2}{2},
```
where ``\mu = \mu(n) = n_\mathrm{max}-n+1``, and ``n_\mathrm{max} = \left\lfloor \lambda \right\rfloor - 1``.
""" E(model::PoschlTeller; n=0)

@doc raw"""
`ψ(model::PoschlTeller, x; n=0)`

```math
\psi_n(x) = P_\lambda^{\mu}(\mathrm{tanh}(x/x_0)) \sqrt{\mu\frac{\Gamma(\lambda-\mu+1)}{\Gamma(\lambda+\mu+1)}},
```
where ``\mu = \mu(n) = n_\mathrm{max}-n+1``, and ``n_\mathrm{max} = \left\lfloor \lambda \right\rfloor - 1`` and ``P_\lambda^{\mu}`` are the associated Legendre functions.
""" ψ(model::PoschlTeller, x; n=0)

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