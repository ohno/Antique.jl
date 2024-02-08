export HydrogenAtom, V, E, ψ, R, L, Y, P

@doc raw"""
HydrogenAtom(Z=1, ℏ=1.0, Eₕ=1.0, a₀=1.0, mₑ=1.0)

``Z`` is the atomic number, ``\hbar`` is the reduced Planck constant (Dirac's constant), ``E_h`` is the Hartree energy, ``a_0``is the Bohr radius and ``m_e`` is the mass of the electron.

""" @kwdef struct HydrogenAtom
  Z = 1
  ℏ = 1.0
  Eₕ = 1.0
  a₀ = 1.0
  mₑ = 1.0
end

@doc raw"""
    V(model::HydrogenAtom, x)

```math
  V(r) = -\frac{Z}{r/a_0}.
```
""" function V(model::HydrogenAtom, r)
  # if r<0
  #   throw(DomainError(r, "r=$r is out of the domain (0≦r)"))
  # end
  Z = model.Z
  a₀ = model.a₀
  return -Z/abs(r/a₀)
end

@doc raw"""
    E(model::HydrogenAtom)

```math
  E = -\frac{Z^2}{2n^2} E_h
```
""" 
function E(model::HydrogenAtom, r, θ, φ; n=1, l=0, m=0)
  Z = model.Z
  Eₕ = model.Eₕ
  return -Z^2/(2*n^2) * Eₕ
end

@doc raw"""
`ψ(model::HydrogenAtom, r, θ, φ; n=1, l=0, m=0)`

```math
    \psi_{nlm}(r) = R_{nl}(r)Y_{lm}(\theta, \phi)
```
"""
function ψ(model::HydrogenAtom, r, θ, φ; n=1, l=0, m=0)
  # if r<0
  #   throw(DomainError(r, "r=$r is out of the domain (0≦r)"))
  # end
  Z = model.Z
  a₀ = model.a₀
  return R(model, r, n=n, l=l) * Y(model, θ, φ, l=l, m=m)
end

@doc raw"""
`R(model::HydrogenAtom, r; n=1, l=0)`

```math
    R_{nl}(r) = -\sqrt{\frac{(n-l-1)!}{2n(n+l)!} \left(\frac{2Z}{n a_0}\right)^3} \left(\frac{2Zr}{n a_0}\right)^l \exp \left(-\frac{Zr}{n a_0}\right) L_{n+l}^{2l+1} \left(\frac{2Zr}{n a_0}\right)
```
""" 
function R(model::HydrogenAtom, r; n=1, l=0)
  # if r<0
  #   throw(DomainError(r, "r=$r is out of the domain (0≦r)"))
  # end
  Z = model.Z
  a₀ = model.a₀
  ρ = 2*Z*abs(r)/(n*a₀)
  N = -sqrt( factorial(n-l-1)/(2*n*factorial(n+l)) * (2*Z/(n*a₀))^3 )
  return N*ρ^l * exp(-ρ/2) * L(model, ρ, n=n+l, k=2*l+1)
end

@doc raw"""
`L(model::HydrogenAtom, x; n=0, k=0)`

```math
    L_{n}^{k}(x) = \sum_{m=0}^{n-k} (-1)^{m+k} \frac{n!}{m!(m+k)!(n-m-k)!} x^m 
```
""" 
function L(model::HydrogenAtom, x; n=0, k=0)
  return sum(m -> (-1)^(m+k) * factorial(n) // (factorial(m) * factorial(m+k) * factorial(n-m-k)) * x^m, 0:n-k)
end      

@doc raw"""
`Y(model::HydrogenAtom, θ, φ; l=0, m=0)`

```math
Y_{lm}(\theta,\varphi) = (-1)^{\frac{|m|+m}{2}} \sqrt{\frac{2l+1}{4\pi} \frac{(l-|m|)!}{(l+|m|)!}} P_l^{|m|} (\cos\theta) \mathrm{e}^{im\varphi}
```
""" 
function Y(model::HydrogenAtom, θ, φ; l=0, m=0)
  N = (im)^(m+abs(m)) * sqrt( (2*l+1)*factorial(l-Int(abs(m))) / (2*factorial(l+Int(abs(m)))) )
  return N * P(model,cos(θ), n=l, m=Int(abs(m))) * exp(im*m*φ) / sqrt(2*π)
end

@doc raw"""
`P(model::HydrogenAtom, x; n=0, m=0)`

```math
\frac{1}{2^n} (1-x^2)^{m/2} \sum_{j=0}^{\left\lfloor\frac{n-m}{2}\right\rfloor} (-1)^j \frac{(2n-2j)!}{j! (n-j)! (n-2j-m)!} x^{(n-2j-m)}
```
""" 
function P(model::HydrogenAtom, x; n=0, m=0)
  return (1//2)^n * (1-x^2)^(m//2) * sum(j -> (-1)^j * factorial(2*n-2*j) // (factorial(j) * factorial(n-j) * factorial(n-2*j-m)) * x^(n-2*j-m), 0:Int(floor((n-m)/2)))
end