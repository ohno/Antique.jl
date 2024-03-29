export PoschlTeller, V, E, ψ

@kwdef struct PoschlTeller
  λ = 1.0
  m = 1.0
  ℏ = 1.0
  x0 = 1.0
end

function V(model::PoschlTeller, x)
  λ = model.λ
  return -λ*(λ+1)/2/cosh(x)^2
end

function nmax(model::PoschlTeller)
  λ = model.λ
  return Int(floor(λ-1)) # if counting n from zero
end

function E(model::PoschlTeller; n=0)
  λ = model.λ
  m = model.m
  ℏ = model.ℏ
  x0 = model.x0

  n_max = nmax(model)
  mu = n_max-n+1
  #m = model.m
  #ℏ = model.ℏ
  if (λ ≈ round(λ)) == false
    @show(λ,round(λ),λ ≈ round(λ))
    error("Error: Currently only integer values for λ are supported.")
  end
  
  if (n ≈ round(n)) == false
    @show(n,round(n),n ≈ round(n))
    error("Error: n is the index of the n-th excited state. It must be an integer.")
  end
  
  if 0 <= n <= n_max
    return -(mu)^2/2 * ℏ^2/(m*x0^2)
  else
    error("Error: n must be non-negative and smaller than λ: 0 <= n < λ")
  end
end

function ψ(model::PoschlTeller, x; n=0)
  λ = model.λ
  n_max = nmax(model)
  mu = n_max-n+1

  if 0 <= n <= n_max
    if λ ≈ round(λ)
      return P(model,tanh(x),n=Int64(λ),m=mu)*sqrt(mu*gamma(λ-mu+1)/gamma(λ+mu+1))
    else
      #throw(DomainError(λ, "Currently only integer values for λ are supported."))
      error("Error: Currently only integer values for λ are supported.")
    end
  else
    #throw(DomainError(n, "n=$n must be non-negative and smaller than λ: 0 <= n < λ"))
    error("Error: n must be non-negative and smaller than λ: 0 <= n < λ")
  end
end

function P(model::PoschlTeller, x; n=0, m=0) # different definition from hydrogen atom: additional factor (-1)^m here
  return (-1)^m * (1//2)^n * (1-x^2)^(m//2) * sum(j -> (-1)^j * factorial(2*n-2*j) // (factorial(j) * factorial(n-j) * factorial(n-2*j-m)) * x^(n-2*j-m), 0:Int(floor((n-m)/2)))
end

# docstrings:
@doc raw"""
`PoschlTeller(λ=1.0,m=1.0,ℏ=1.0,x0=1.0)`

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
`nmax(model::PoschlTeller)`

```math
n_\mathrm{max} = \left\lfloor \lambda \right\rfloor - 1.
```
""" nmax(model::PoschlTeller)

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

Associated Legendre polynomials are the associated Legendre functions for integer indices. Please note here, that for the Poschl-Teller potential we use a slightly different notation of the associated Legendre functions as compared to the model HydrogenAtom. Here we have an additional factor ``(-1)^\m``.

""" P(model::PoschlTeller, x; n=0, m=0)