export InfinitePotentialWell3D, V, E, ψ

# parameters
@kwdef struct InfinitePotentialWell3D
  L = [1.0, 1.0, 1.0]
  m = 1.0
  ℏ = 1.0
end

# potential
function V(model::InfinitePotentialWell3D, x)
  L = model.L
  return prod(@. 0<x<L) ? 0 : Inf
end

# eigenvalues
function E(model::InfinitePotentialWell3D; n::Vector{Int}=[1,1,1])
  if !(prod(1 .≤ n))
    throw(DomainError("n = $n", "This function is defined for 1 .≤ n."))
  end
  L = model.L
  m = model.m
  ℏ = model.ℏ
  return sum(@. (ℏ^2*n^2*π^2) / (2*m*L^2))
end

# eigenfunctions
function ψ(model::InfinitePotentialWell3D, x; n::Vector{Int}=[1,1,1])
  if !(prod(1 .≤ n))
    throw(DomainError("n = $n", "This function is defined for 1 .≤ n."))
  end
  L = model.L
  return prod(@. 0<x<L) ? prod(@. sqrt(2/L) * sin(n*π*x/L)) : 0
end

# docstrings

@doc raw"""
## Model

This model is described with the time-independent Schrödinger equation
```math
  \hat{H} \psi(x,y,z) = E \psi(x,y,z),
```
and the Hamiltonian
```math
  \hat{H} = - \frac{\hbar^2}{2m} \left(\frac{\partial^2}{\partial x ^2} + \frac{\partial^2}{\partial y ^2} + \frac{\partial^2}{\partial z ^2}\right) + V(x,y,z).
```
Parameters are specified with the following struct:

```
IPW3D = InfinitePotentialWell3D(L=[1.0,1.0,1.0], m=1.0, ℏ=1.0)
```

``L`` is a vector of the lengths of the box in ``x``,``y``,``z``-direction, ``m`` is the mass of the particle and ``\hbar`` is the reduced Planck constant (Dirac's constant).

## References

* [McQuarrie1997](@cite) D. A. McQuarrie, J. D. Simon, _Physical chemistry : a molecular approach_ (University Science Books, 1997), (https://mitpress.mit.edu/9781940380216/physical-chemistry/) p.90, 3-9. The Problem of a Particle in a Three-Dimensional Box Is a Simple Extension of the One-Dimensional Case
""" InfinitePotentialWell3D

@doc raw"""
`V(model::InfinitePotentialWell3D, x)`

```math
V(x,y,z) =
\left\{
  \begin{array}{ll}
  0      & 0 \leq x \leq L_x \ \mathrm{and}\  0 \leq y \leq L_y \ \mathrm{and}\  0 \leq z \leq L_z \\
  \infty & \mathrm{elsewhere}
  \end{array}
\right.
```
""" V(model::InfinitePotentialWell3D, x)

@doc raw"""
`E(model::InfinitePotentialWell3D; n::Vector{Int}=[1,1,1])`

```math
E_{n_x,n_y,n_z}
= \frac{\hbar^2 n_x^2 \pi^2}{2 m L_x^2}
+ \frac{\hbar^2 n_y^2 \pi^2}{2 m L_y^2}
+ \frac{\hbar^2 n_z^2 \pi^2}{2 m L_z^2}
```
""" E(model::InfinitePotentialWell3D; n::Vector{Int})

@doc raw"""
`ψ(model::InfinitePotentialWell3D, x; n::Vector{Int}=[1,1,1])`

The wave functions can be expressed as products of wave functions in a one-dimensional box.

```math
\begin{aligned}
  \psi_{n_x,n_y,n_z}(x,y,z)
  &=     \psi_{n_x}(x)
  \times \psi_{n_y}(y)
  \times \psi_{n_z}(z) \\
  &=     \sqrt{\frac{2}{L_x}} \sin \frac{n_x \pi x}{L_x}
  \times \sqrt{\frac{2}{L_y}} \sin \frac{n_y \pi y}{L_y}
  \times \sqrt{\frac{2}{L_z}} \sin \frac{n_z \pi z}{L_z}
\end{aligned}
```
""" ψ(model::InfinitePotentialWell3D, x; n::Vector{Int})