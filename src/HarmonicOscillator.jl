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
`HarmonicOscillator(k=1.0, m=1.0, ℏ=1.0)`

``k`` is the force constant, ``m`` is the mass of particle and ``\hbar`` is the reduced Planck constant (Dirac's constant).
""" HarmonicOscillator

@doc raw"""
`V(model::HarmonicOscillator, x)`

```math
V(x)
= \frac{1}{2} k x^2
= \frac{1}{2} m \omega^2 x^2
= \frac{1}{2} \hbar \omega \xi^2,
```
where ``\omega = \sqrt{k/m}`` is the angular frequency and ``\xi = \sqrt{\frac{m\omega}{\hbar}}x``.
""" V(model::HarmonicOscillator, x)

@doc raw"""
`E(model::HarmonicOscillator; n=0)`

```math
E_n = \hbar \omega \left( n + \frac{1}{2} \right),
```
where ``\omega = \sqrt{k/m}`` is the angular frequency.
""" E(model::HarmonicOscillator; n=0)

@doc raw"""
`ψ(model::HarmonicOscillator, x; n=0)`

```math
\psi_n(x) = A_n H_n(\xi) \exp{\left( -\frac{\xi^2}{2} \right)},
```
where ``\omega = \sqrt{k/m}``, ``\xi = \sqrt{\frac{m\omega}{\hbar}}x``, ``A_n = \sqrt{\frac{1}{n! 2^n} \sqrt{\frac{m\omega}{\pi\hbar}}}``, ``H_n(x) = (-1)^n \mathrm{e}^{x^2} \frac{\mathrm{d}^n}{\mathrm{d}x^n} \mathrm{e}^{-x^2}`` are defined.
""" ψ(model::HarmonicOscillator, x; n=0)

@doc raw"""
`H(model::HarmonicOscillator, x; n=0)`

Rodrigues' formula & closed-form:
```math
\begin{aligned}
  H_{n}(x)
  &:= (-1)^n \mathrm{e}^{x^2} \frac{\mathrm{d}^n}{\mathrm{d}x^n} \mathrm{e}^{-x^2} \\
  &= n! \sum_{m=0}^{\lfloor n/2 \rfloor} \frac{(-1)^m}{m! (n-2m)!}(2 x)^{n-2m}.
\end{aligned}
```
Examples:
```math
\begin{aligned}
  H_{0}(x)  &= 1, \\
  H_{1}(x)  &= 2 x, \\
  H_{2}(x)  &= -2 + 4 x^{2}, \\
  H_{3}(x)  &= -12 x + 8 x^{3}, \\
  H_{4}(x)  &= 12 - 48 x^{2} + 16 x^{4}, \\
  H_{5}(x)  &= 120 x - 160 x^{3} + 32 x^{5}, \\
  H_{6}(x)  &= -120 + 720 x^{2} - 480 x^{4} + 64 x^{6}, \\
  H_{7}(x)  &= -1680 x + 3360 x^{3} - 1344 x^{5} + 128 x^{7}, \\
  H_{8}(x)  &= 1680 - 13440 x^{2} + 13440 x^{4} - 3584 x^{6} + 256 x^{8}, \\
  H_{9}(x)  &= 30240 x - 80640 x^{3} + 48384 x^{5} - 9216 x^{7} + 512 x^{9}, \\
  &\vdots
\end{aligned}
```
""" H(model::HarmonicOscillator, x; n=0)
