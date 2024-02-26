export PoschlTeller, V, E, ψ, H

# PT = PoschlTeller(lambda=2.0) does not work. "UndefVarError: PoschlTeller not defined"

@kwdef struct PoschlTeller
  lambda = 1.0
  #m = 1.0
  #ℏ = 1.0
end

function V(model::PoschlTeller, x)
  lambda = model.lambda
  return -lambda*(lambda+1)/2/cosh(x)^2
end

function nmax(model::PoschlTeller)
  lambda = model.lambda
  return Int(floor(lambda-1)) # if counting n from zero
end

function E(model::PoschlTeller; n=0)
  lambda = model.lambda
  n_max = nmax(model)
  mu = n_max-n+1
  #m = model.m
  #ℏ = model.ℏ
  if (lambda ≈ round(lambda)) == false
    @show(lambda,round(lambda),lambda ≈ round(lambda))
    error("Error: Currently only integer values for lambda are supported.")
  end
  
  if (n ≈ round(n)) == false
    @show(n,round(n),n ≈ round(n))
    error("Error: n is the index of the n-th excited state. It must be an integer.")
  end
  
  if 0 <= n <= n_max
    return -(mu)^2/2
  else
    error("Error: n must be non-negative and smaller than lambda: 0 <= n < lambda")
  end
end

function ψ(model::PoschlTeller, x; n=0)
  lambda = model.lambda
  n_max = nmax(model)
  mu = n_max-n+1

  if 0 <= n <= n_max
    if lambda ≈ round(lambda)
      return P(model,tanh(x),n=Int64(lambda),m=mu)*sqrt(mu*gamma(lambda-mu+1)/gamma(lambda+mu+1))
    else
      #throw(DomainError(lambda, "Currently only integer values for lambda are supported."))
      error("Error: Currently only integer values for lambda are supported.")
    end
  else
    #throw(DomainError(n, "n=$n must be non-negative and smaller than lambda: 0 <= n < lambda"))
    error("Error: n must be non-negative and smaller than lambda: 0 <= n < lambda")
  end
end

function P(model::PoschlTeller, x; n=0, m=0) # different definition from hydrogen atom: additional factor (-1)^m here
  return (-1)^m * (1//2)^n * (1-x^2)^(m//2) * sum(j -> (-1)^j * factorial(2*n-2*j) // (factorial(j) * factorial(n-j) * factorial(n-2*j-m)) * x^(n-2*j-m), 0:Int(floor((n-m)/2)))
end

# docstrings:
@doc raw"""
`PoschlTeller(lambda=1.0)`

``\lambda`` determines the potential strength. This model is defined dimensionless, i.e. ``x`` is given in units of a characteristic length ``x_0``, and ``E`` in units of a characteristic energy, e.g. ``E_\mathrm{char} = \frac{\hbar^2}{2 m x_0^2}``.
""" PoschlTeller

@doc raw"""
`V(model::PoschlTeller, x)`

```math
\begin{aligned}
  V(x)
  &= -\frac{\lambda(\lambda+1)}{2}  \mathrm{sech}(x)^2
  &= -\frac{\lambda(\lambda+1)}{2}  \frac{1}{\mathrm{cosh}(x)^2}.
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
E_n = -\frac{\mu^2}{2},
```
where ``\mu = \mu(n) = n_\mathrm{max}-n+1``, and ``n_\mathrm{max} = \left\lfloor \lambda \right\rfloor - 1``.
""" E(model::PoschlTeller; n=0)

@doc raw"""
`ψ(model::PoschlTeller, x; n=0)`

```math
\psi_n(x) = P_\lambda^{\mu}(\mathrm{tanh}(x)) \sqrt{\mu\frac{\Gamma(\lambda-\mu+1)}{\Gamma(\lambda+\mu+1)}},
```
where ``\mu = \mu(n) = n_\mathrm{max}-n+1``, and ``n_\mathrm{max} = \left\lfloor \lambda \right\rfloor - 1`` and ``P_\lambda^{\mu}`` are the associated Legendre functions.
""" ψ(model::PoschlTeller, x; n=0)

@doc raw"""
`P(model::PoschlTeller, x; n=0, m=0)`

Associated Legendre polynomials are the associated Legendre functions for integer indices. Please note here, that for the Poschl-Teller potential we use a slightly different notation of the associated Legendre functions as compared to the model HydrogenAtom. Here we have an additional factor ``(-1)^\m``.

""" P(model::PoschlTeller, x; n=0, m=0)