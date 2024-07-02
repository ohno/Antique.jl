export InfinitePotentialWell, V, E, ψ

# parameters
@kwdef struct InfinitePotentialWell
  L = 1.0
  m = 1.0
  ℏ = 1.0
end

# potential
function V(model::InfinitePotentialWell, x)
  L = model.L
  return 0<x<L ? 0 : Inf
end

# eigenvalues
function E(model::InfinitePotentialWell; n::Int=1)
  if !(1 ≤ n)
    throw(DomainError("n = $n", "n must be 1 or more: 1 ≤ n."))
  end
  L = model.L
  m = model.m
  ℏ = model.ℏ
  return (ℏ^2*n^2*π^2) / (2*m*L^2)
end

# eigenfunctions
function ψ(model::InfinitePotentialWell, x; n::Int=1)
  if !(1 ≤ n)
    throw(DomainError("n = $n", "n must be 1 or more: 1 ≤ n."))
  end
  L = model.L
  return 0<x<L ? sqrt(2/L) * sin(n*π*x/L) : 0
end

# docstrings

@doc raw"""
This model is described with the time-independent Schrödinger equation
```math
  \hat{H} \psi(x) = E \psi(x),
```
and the Hamiltonian
```math
  \hat{H} = - \frac{\hbar^2}{2m} \frac{\mathrm{d}^2}{\mathrm{d}x ^2} + V(x).
```
Parameters are specified with the following struct:

`InfinitePotentialWell(L=1.0, m=1.0, ℏ=1.0)`

``L`` is the length of the box, ``m`` is the mass of particle and ``\hbar`` is the reduced Planck constant (Dirac's constant).

Reference:
- [D. J. Griffiths, D. F. Schroeter, _Introduction to Quantum Mechanics_ **Third Edition** (Cambridge University Press, 2018)](https://doi.org/10.1017/9781316995433) p.31, 2.2 THE INFINITE SQUARE WELL

Proofs:
- [Eigen Functions & Eigen Values](https://ja.wolframalpha.com/input?i2d=true&i=D%5B%5C%2840%29Sqrt%5BDivide%5B2%2Ca%5D%5Dsin%5C%2840%29Divide%5Bn%CF%80x%2Ca%5D%5C%2841%29%5C%2841%29%2C%7Bx%2C2%7D%5D)
- [Normalization](https://ja.wolframalpha.com/input?i=Integrate%5B%28%28Sqrt%5B2%2Fa%5Dsin%28%CF%80x%2Fa%29%29%29%5E2%2C+%7Bx%2C0%2Ca%7D%5D)
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
`E(model::InfinitePotentialWell; n::Int=1)`

```math
E_n = \frac{\hbar^2 n^2 \pi^2}{2 m L^2}
```
""" E(model::InfinitePotentialWell; n::Int=1)

@doc raw"""
`ψ(model::InfinitePotentialWell, x; n::Int=1)`

```math
\psi_n(x) = \sqrt{\frac{2}{L}} \sin \frac{n\pi x}{L}
```
""" ψ(model::InfinitePotentialWell, x; n::Int=1)