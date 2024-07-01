export CoulombTwoBody, V, E, ψ, R, L, Y, P

# parameters
@kwdef struct CoulombTwoBody
  z₁ = -1
  z₂ = 1
  m₁ = 1.0
  m₂ = 1.0
  mₑ = 1.0
  a₀ = 1.0
  Eₕ = 1.0
  ℏ = 1.0
end

# potential
function V(model::CoulombTwoBody, r)
  if !(0 ≤ r)
    throw(DomainError("r = $r", "r must be non-negative: 0 ≤ r."))
  end
  z₁ = model.z₁
  z₂ = model.z₂
  a₀ = model.a₀
  Eₕ = model.Eₕ
  return z₁*z₂/abs(r/a₀) * Eₕ
end

# eigenvalues
function E(model::CoulombTwoBody; n::Int=1)
  if !(1 ≤ n)
    throw(DomainError("n = $n", "n must be 1 or more: 1 ≤ n."))
  end
  z₁ = model.z₁
  z₂ = model.z₂
  m₁ = model.m₁
  m₂ = model.m₂
  mₑ = model.mₑ
  μ = (1/m₁ + 1/m₂)^(-1)
  Eₕ = model.Eₕ
  return -(z₁*z₂)^2/(2*n^2) * μ/mₑ * Eₕ
end

# eigenfunctions
function ψ(model::CoulombTwoBody, r, θ, φ; n::Int=1, l::Int=0, m::Int=0)
  if !(1 ≤ n && 0 ≤ l < n && -l ≤ m ≤ l)
    throw(DomainError("(n,l,m) = ($n,$l,$m)", "This function is defined for 1 ≤ n, 0 ≤ l < n and -l ≤ m ≤ l."))
  end
  if !(0 ≤ r && 0 ≤ θ < π && 0 ≤ φ < 2π)
    throw(DomainError("(r,θ,φ) = ($r,$θ,$φ)", "This function is defined for 0 ≤ r, 0 ≤ θ < π, 0 ≤ φ < 2π."))
  end
  return R(model, r, n=n, l=l) * Y(model, θ, φ, l=l, m=m)
end

# radial function
function R(model::CoulombTwoBody, r; n=1, l=0)
  z₁ = model.z₁
  z₂ = model.z₂
  a₀ = model.a₀
  m₁ = model.m₁
  m₂ = model.m₂
  mₑ = model.mₑ
  μ  = (1/m₁ + 1/m₂)^(-1)
  aμ = a₀ * mₑ / μ
  ρ = 2*(-z₁*z₂)*abs(r)/(n*aμ)
  N = -sqrt( factorial(n-l-1)/(2*n*factorial(n+l)) * (2*(-z₁*z₂)/(n*aμ))^3 )
  return N*ρ^l * exp(-ρ/2) * L(model, ρ, n=n+l, k=2*l+1)
end

# associated Laguerre polynomials
function L(model::CoulombTwoBody, x; n=0, k=0)
  return sum(m -> (-1)^(m+k) * factorial(n) // (factorial(m) * factorial(m+k) * factorial(n-m-k)) * x^m, 0:n-k)
end      

# spherical harmonics
function Y(model::CoulombTwoBody, θ, φ; l=0, m=0)
  N = (-1)^((abs(m)+m)/2) * sqrt( (2*l+1)*factorial(l-Int(abs(m))) / (2*factorial(l+Int(abs(m)))) )
  return N * P(model,cos(θ), n=l, m=Int(abs(m))) * exp(im*m*φ) / sqrt(2*π)
end

# associated Legendre polynomials
function P(model::CoulombTwoBody, x; n=0, m=0)
  return (1//2)^n * (1-x^2)^(m//2) * sum(j -> (-1)^j * factorial(2*n-2*j) // (factorial(j) * factorial(n-j) * factorial(n-2*j-m)) * x^(n-2*j-m), 0:Int(floor((n-m)/2)))
end

# docstrings

@doc raw"""
`CoulombTwoBody(z₁=-1, z₂=1, m₁=1.0, m₂=1.0, mₑ=1.0, a₀=1.0, Eₕ=1.0, ℏ=1.0)`

``z₁`` is the charge number of particle 1, 
``z₂`` is the charge number of particle 2, 
``m₁`` is the mass of particle 1, 
``m₂`` is the mass of particle 2,
``m_\mathrm{e}`` is the electron mass (use the same unit as ``m₁`` and ``m₂``. For example of hydrogen atom, use ``m_\mathrm{e}=9.1093837139\times10^{-31}\mathrm{kg}``, ``m_1=9.1093837139\times10^{-31}\mathrm{kg}`` and ``m_2=1.67262192595\times10^{-27}\mathrm{kg}`` in the IS unit system, use ``~m_\mathrm{e}=1.0~m_\mathrm{e}``, ``m_1=1.0~m_\mathrm{e}`` and ``m_2=1836.152673426~m_\mathrm{e}`` in the atomic unit.),
``a_0`` is the Bohr radius,
``E_\mathrm{h}`` is the Hartree energy and
``\hbar`` is the reduced Planck constant (Dirac's constant).
""" CoulombTwoBody

@doc raw"""
`V(model::CoulombTwoBody, r)`

```math
\begin{aligned}
  V(r)
  &= - \frac{Ze^2}{4\pi\varepsilon_0 r} 
  &= - \frac{e^2}{4\pi\varepsilon_0 a_0} \frac{Z}{r/a_0}
  &= - \frac{Z}{r/a_0} E_\mathrm{h},
\end{aligned}
```
where ``E_\mathrm{h} = \frac{\hbar^2}{m_\mathrm{e}{a_0}^2} = \frac{e^2}{4\pi\varepsilon_0a_0} = \frac{m_\mathrm{e}e^4}{\left(4\pi\varepsilon_0\right)^2\hbar^2}`` is the Hartree energy, one of atomic unit. The domain is $0\leq r \lt \infty$.
""" V(model::CoulombTwoBody, r)

@doc raw"""
`E(model::CoulombTwoBody; n::Int=1)`

```math
E_n
= -\frac{(z_1 z_2)^2}{2n^2} \frac{\mu}{m_\mathrm{e}} E_\mathrm{h},
```
where $\mu=\left(\frac{1}{m_1}+\frac{1}{m_2}\right)^{-1}$ is the reduced mass of particle 1 and particle 2, ``E_\mathrm{h} = \frac{\hbar^2}{m_\mathrm{e}{a_0}^2} = \frac{e^2}{4\pi\varepsilon_0a_0} = \frac{m_\mathrm{e}e^4}{\left(4\pi\varepsilon_0\right)^2\hbar^2}`` is the Hartree energy, one of atomic unit. About atomic units, see section 3.9.2 of the [IUPAC GreenBook](https://iupac.org/what-we-do/books/greenbook/). In other units, ``E_\mathrm{h} = 27.211~386~245~988(53)~\mathrm{eV}`` from [here](https://physics.nist.gov/cgi-bin/cuu/Value?hrev).
""" E(model::CoulombTwoBody; n::Int=1)

@doc raw"""
`ψ(model::CoulombTwoBody, r, θ, φ; n::Int=1, l::Int=0, m::Int=0)`

```math
\psi_{nlm}(\pmb{r}) = R_{nl}(r) Y_{lm}(\theta,\varphi)
```
The domain is $0\leq r \lt \infty, 0\leq \theta \lt \pi, 0\leq \varphi \lt 2\pi$.
""" ψ(model::CoulombTwoBody, r, θ, φ; n::Int=1, l::Int=0, m::Int=0)

@doc raw"""
`R(model::CoulombTwoBody, r; n=1, l=0)`

```math
R_{nl}(r) = -\sqrt{\frac{(n-l-1)!}{2n(n+l)!} \left(\frac{2Z}{n a_\mu}\right)^3} \left(\frac{2Zr}{n a_\mu}\right)^l \exp \left(-\frac{Zr}{n a_\mu}\right) L_{n+l}^{2l+1} \left(\frac{2Zr}{n a_\mu}\right),
```
where ``\frac{1}{\mu} = \frac{1}{m_1}+\frac{1}{m_2}``, ``a_\mu = a_0 \frac{m_\mathrm{e}}{\mu}``, Laguerre polynomials are defined as ``L_n(x) = \frac{1}{n!} \mathrm{e}^x \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left( \mathrm{e}^{-x} x^n \right)``, and associated Laguerre polynomials are defined as ``L_n^{k}(x) = \frac{\mathrm{d}^k}{\mathrm{d}x^k} L_n(x)``. Note that replace ``2n(n+l)!`` with ``2n[(n+l)!]^3`` if Laguerre polynomials are defined as ``L_n(x) = \mathrm{e}^x \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left( \mathrm{e}^{-x} x^n \right)``. 
The domain is $0\leq r \lt \infty$.
""" R(model::CoulombTwoBody, r; n=1, l=0)

@doc raw"""
`L(model::CoulombTwoBody, x; n=0, k=0)`

Rodrigues' formula & closed-form:
```math
\begin{aligned}
L_n^{k}(x)
  &= \frac{\mathrm{d}^k}{\mathrm{d}x^k} L_n(x) \\
  &= \frac{\mathrm{d}^k}{\mathrm{d}x^k} \frac{1}{n!} \mathrm{e}^x \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left( \mathrm{e}^{-x} x^n \right) \\
  &= \sum_{m=0}^{n-k} (-1)^{m+k} \frac{n!}{m!(m+k)!(n-m-k)!} x^m \\
  &= (-1)^k L_{n-k}^{(k)}(x),
\end{aligned}
```
where Laguerre polynomials are defined as ``L_n(x)=\frac{1}{n!}\mathrm{e}^x \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left( \mathrm{e}^{-x} x^n \right)``.

Examples:
```math
\begin{aligned}
  L_0^0(x) &= 1, \\
  L_1^0(x) &= 1 - x, \\
  L_1^1(x) &= 1, \\
  L_2^0(x) &= 1 - 2 x + 1/2 x^2, \\
  L_2^1(x) &= 2 - x, \\
  L_2^2(x) &= 1, \\
  L_3^0(x) &= 1 - 3 x + 3/2 x^2 - 1/6 x^3, \\
  L_3^1(x) &= 3 - 3 x + 1/2 x^2, \\
  L_3^2(x) &= 3 - x, \\
  L_3^3(x) &= 1, \\
  L_4^0(x) &= 1 - 4 x + 3 x^2 - 2/3 x^3 + 5/12 x^4, \\
  L_4^1(x) &= 4 - 6 x + 2 x^2 - 1/6 x^3, \\
  L_4^2(x) &= 6 - 4 x + 1/2 x^2, \\
  L_4^3(x) &= 4 - x, \\
  L_4^4(x) &= 1, \\
  \vdots
\end{aligned}
```
""" L(model::CoulombTwoBody, x; n=0, k=0)

@doc raw"""
`Y(model::CoulombTwoBody, θ, φ; l=0, m=0)`

```math
Y_{lm}(\theta,\varphi) = (-1)^{\frac{|m|+m}{2}} \sqrt{\frac{2l+1}{4\pi} \frac{(l-|m|)!}{(l+|m|)!}} P_l^{|m|} (\cos\theta) \mathrm{e}^{im\varphi}.
```
The domain is $0\leq \theta \lt \pi, 0\leq \varphi \lt 2\pi$. Note that some variants are connected by 
```math
i^{|m|+m} \sqrt{\frac{(l-|m|)!}{(l+|m|)!}} P_l^{|m|} = (-1)^{\frac{|m|+m}{2}} \sqrt{\frac{(l-|m|)!}{(l+|m|)!}} P_l^{|m|} = (-1)^m \sqrt{\frac{(l-m)!}{(l+m)!}} P_l^{m}.
```

""" Y(model::CoulombTwoBody, θ, φ; l=0, m=0)

@doc raw"""
`P(model::CoulombTwoBody, x; n=0, m=0)`

Rodrigues' formula & closed-form:
```math
\begin{aligned}
  P_n^m(x)
  &= \left( 1-x^2 \right)^{m/2} \frac{\mathrm{d}^m}{\mathrm{d}x^m} P_n(x) \\
  &= \left( 1-x^2 \right)^{m/2} \frac{\mathrm{d}^m}{\mathrm{d}x^m} \frac{1}{2^n n!} \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left[ \left( x^2-1 \right)^n \right] \\
  &= \frac{1}{2^n} (1-x^2)^{m/2} \sum_{j=0}^{\left\lfloor\frac{n-m}{2}\right\rfloor} (-1)^j \frac{(2n-2j)!}{j! (n-j)! (n-2j-m)!} x^{(n-2j-m)}.
\end{aligned},
```
where Legendre polynomials are defined as ``P_n(x) = \frac{1}{2^n n!} \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left[ \left( x^2-1 \right)^n \right]``. Note that ``P_l^{-m} = (-1)^m \frac{(l-m)!}{(l+m)!} P_l^m`` for ``m<0``. (It is not compatible with ``P_k^m(t) = (-1)^m\left( 1-t^2 \right)^{m/2} \frac{\mathrm{d}^m P_k(t)}{\mathrm{d}t^m}`` caused by ``(-1)^m``.) The specific formulae are given below.

Examples:
```math
\begin{aligned}
  P_{0}^{0}(x) &= 1, \\
  P_{1}^{0}(x) &= x, \\
  P_{1}^{1}(x) &= \left(+1\right)\sqrt{1-x^2}, \\
  P_{2}^{0}(x) &= -1/2 + 3/2 x^{2}, \\
  P_{2}^{1}(x) &= \left(-3 x\right)\sqrt{1-x^2}, \\
  P_{2}^{2}(x) &= 3 - 6 x, \\
  P_{3}^{0}(x) &= -3/2 x + 5/2 x^{3}, \\
  P_{3}^{1}(x) &= \left(3/2 - 15/2 x^{2}\right)\sqrt{1-x^2}, \\
  P_{3}^{2}(x) &= 15 x - 30 x^{2}, \\
  P_{3}^{3}(x) &= \left(15 - 30 x\right)\sqrt{1-x^2}, \\
  P_{4}^{0}(x) &= 3/8 - 15/4 x^{2} + 35/8 x^{4}, \\
  P_{4}^{1}(x) &= \left(- 15/2 x + 35/2 x^{3}\right)\sqrt{1-x^2}, \\
  P_{4}^{2}(x) &= -15/2 + 15 x + 105/2 x^{2} - 105 x^{3}, \\
  P_{4}^{3}(x) &= \left(105 x - 210 x^{2}\right)\sqrt{1-x^2}, \\
  P_{4}^{4}(x) &= 105 - 420 x + 420 x^{2}, \\
  & \vdots
\end{aligned}
```
""" P(model::CoulombTwoBody, x; n=0, m=0)
