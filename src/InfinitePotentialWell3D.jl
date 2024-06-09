export InfinitePotentialWell3D, V, E, ψ

# parameters
@kwdef struct InfinitePotentialWell3D
  Lx = 1.0
  Ly = 1.0
  Lz = 1.0
  m = 1.0
  ℏ = 1.0
end

# potential
function V(model::InfinitePotentialWell3D, x,y,z)
  Lx = model.Lx
  Ly = model.Ly
  Lz = model.Lz
  return (0<x<Lx&&0<y<Ly&&0<z<Lz) ? 0 : Inf
end

# eigenvalues
function E(model::InfinitePotentialWell3D; nx::Int=1, ny::Int=1, nz::Int=1)
  if !(1 ≤ nx && 1 ≤ ny && 1 ≤ nz)
    throw(DomainError("(nx,ny,nz) = ($nx,$ny,$nz)", "This function is defined for 1 ≤ nx, 1 ≤ ny and 1 ≤ nz."))
  end
  Lx = model.Lx
  Ly = model.Ly
  Lz = model.Lz
  m = model.m
  ℏ = model.ℏ
  return (ℏ^2*π^2) / (2*m) * (nx^2/Lx^2 + ny^2/Ly^2 + nz^2/Lz^2)
end

# eigenfunctions
function ψ(model::InfinitePotentialWell3D, x,y,z; nx::Int=1, ny::Int=1, nz::Int=1)
  if !(1 ≤ nx && 1 ≤ ny && 1 ≤ nz)
    throw(DomainError("(nx,ny,nz) = ($nx,$ny,$nz)", "This function is defined for 1 ≤ nx, 1 ≤ ny and 1 ≤ nz."))
  end
  Lx = model.Lx
  Ly = model.Ly
  Lz = model.Lz
  return (0<x<Lx&&0<y<Ly&&0<z<Lz) ? sqrt(8/(Lx*Ly*Lz))*sin(nx*π*x/Lx)*sin(ny*π*y/Ly)*sin(nz*π*z/Lz) : 0
end

# docstrings

@doc raw"""
`InfinitePotentialWell3D(Lx=1.0, Ly=1.0, Lz=1.0, m=1.0, ℏ=1.0)`

``L_x,L_y,L_z`` are the lengths of the box in ``x``,``y``,``z``-direction, ``m`` is the mass of particle and ``\hbar`` is the reduced Planck constant (Dirac's constant).
""" InfinitePotentialWell3D

@doc raw"""
`V(model::InfinitePotentialWell3D, x,y,z)`

```math
V(x,y,z) =
\left\{
  \begin{array}{ll}
  0      & 0 \leq x \leq L_x \ \mathrm{and}\  0 \leq y \leq L_y \ \mathrm{and}\  0 \leq z \leq L_z \\
  \infty & \mathrm{elsewhere}
  \end{array}
\right.
```
""" V(model::InfinitePotentialWell3D, x,y,z)

@doc raw"""
`E(model::InfinitePotentialWell3D; nx=1, ny=1, nz=1)`

```math
E_{n_x,n_y,n_z} = \frac{\hbar^2 \pi^2}{2 m} \left(\frac{n_x^2}{L_x^2} + \frac{n_y^2}{L_y^2} + \frac{n_z^2}{L_z^2}\right)
```
""" E(model::InfinitePotentialWell3D; nx=1, ny=1, nz=1)

@doc raw"""
`ψ(model::InfinitePotentialWell3D, x,y,z; nx=1, ny=1, nz=1)`

The wave functions can be expressed as products of wave functions in a one-dimensional box.

```math
\psi_{n_x,n_y,n_z}(x,y,z) = \psi_{n_x}(x)\psi_{n_y}(y)\psi_{n_z}(z) = \sqrt{\frac{8}{L_xL_yL_z}} \sin\left(\frac{n_x\pi x}{L_x}\right) \sin\left(\frac{n_y\pi y}{L_y}\right) \sin\left(\frac{n_z\pi z}{L_z}\right)
```
""" ψ(model::InfinitePotentialWell3D, x,y,z; nx=1, ny=1, nz=1)