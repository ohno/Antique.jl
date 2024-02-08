export HarmonicOscillator, V, E, ψ, H


@kwdef struct HarmonicOscillator
  k = 1.0
  m = 1.0
  ℏ = 1.0
end

function V(model::HarmonicOscillator, x)
  k = model.k
  return 1//2 * k * x^2
end

function E(model::HarmonicOscillator; n=0)
  k = model.k
  m = model.m
  ℏ = model.ℏ
  ω = sqrt(k/m)
  return ℏ * ω * (n+1//2)
end


function ψ(model::HarmonicOscillator, x; n=0)
  k = model.k
  m = model.m
  ℏ = model.ℏ
  ω = sqrt(k/m)
  A = sqrt(1//(factorial(n)*2^n)*sqrt(m*ω/(π*ℏ)))
  ξ = sqrt(m*ω/ℏ) * x
  return A * H(model,ξ,n=n) * exp(-ξ^2/2)
end


function H(model::HarmonicOscillator, x; n=0)
  return factorial(n) * sum(i -> (-1)^i // (factorial(i)  * factorial(n-2*i)) * (2*x)^(n-2*i), 0:Int(floor(n/2)))
end

@doc raw"""
HarmonicOscillator(α=1.0, m=1.0, ℏ=1.0)

``k`` is the force constant, ``m`` is the mass of particle and ``\hbar`` is the reduced Planck constant (Dirac's constant).
""" HarmonicOscillator

@doc raw"""
    V(model::HarmonicOscillator, x)

```math
  V(x) = \frac{1}{2} k x^2.
```
""" V(model::HarmonicOscillator, x)

@doc raw"""
    E(model::HarmonicOscillator)

```math
  E = -\hbar \omega \left(n+\frac{1}{2}\right)
```
""" E(model::HarmonicOscillator; n=0)

@doc raw"""
`ψ(model::HarmonicOscillator, x; n=0)`

```math
   \psi(x) = (2^n n!)^{-1/2} \left(\frac{m \omega}{\pi \hbar} \right)^{1/4} \mathrm{e}^{\xi^2 /2} H_n(\xi)
```
""" ψ(model::HarmonicOscillator, x; n=0)

@doc raw"""
`H(model::HarmonicOscillator, x; n=0)`

```math
H_n(x) = n! \sum_{m=0}^{n} \frac{(-1)^m}{m!(n - 2m)!} (2x)^{n-2m}
```
""" H(model::HarmonicOscillator, x; n=0)