export HydrogenAtom, V, E, ψ, R, L, Y, P

@kwdef struct HydrogenAtom
  Z = 1
  mₑ = 1.0
  a₀ = 1.0
  Eₕ = 1.0
  ℏ = 1.0
end

function V(model::HydrogenAtom, r)
  # if r<0
  #   throw(DomainError(r, "r=$r is out of the domain (0≦r)"))
  # end
  Z = model.Z
  a₀ = model.a₀
  Eₕ = model.Eₕ
  return Eₕ*-1*Z/abs(r/a₀)
end

function E(model::HydrogenAtom; n=1)
  Z = model.Z
  Eₕ = model.Eₕ
  return -Z^2/(2*n^2) * Eₕ
end

function ψ(model::HydrogenAtom, r, θ, φ; n=1, l=0, m=0)
  # if r<0
  #   throw(DomainError(r, "r=$r is out of the domain (0≦r)"))
  # end
  return R(model, r, n=n, l=l) * Y(model, θ, φ, l=l, m=m)
end

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

function L(model::HydrogenAtom, x; n=0, k=0)
  return sum(m -> (-1)^(m+k) * factorial(n) // (factorial(m) * factorial(m+k) * factorial(n-m-k)) * x^m, 0:n-k)
end      

function Y(model::HydrogenAtom, θ, φ; l=0, m=0)
  N = (im)^(m+abs(m)) * sqrt( (2*l+1)*factorial(l-Int(abs(m))) / (2*factorial(l+Int(abs(m)))) )
  return N * P(model,cos(θ), n=l, m=Int(abs(m))) * exp(im*m*φ) / sqrt(2*π)
end

function P(model::HydrogenAtom, x; n=0, m=0)
  return (1//2)^n * (1-x^2)^(m//2) * sum(j -> (-1)^j * factorial(2*n-2*j) // (factorial(j) * factorial(n-j) * factorial(n-2*j-m)) * x^(n-2*j-m), 0:Int(floor((n-m)/2)))
end

@doc raw"""
`HydrogenAtom(Z=1, mₑ=1.0, a₀=1.0, Eₕ=1.0, ℏ=1.0)`

``Z`` is the atomic number, ``m_\mathrm{e}`` is the electron mass, ``a_0``is the Bohr radius, ``E_\mathrm{h}`` is the Hartree energy and ``\hbar`` is the reduced Planck constant (Dirac's constant).
""" HydrogenAtom

@doc raw"""
`V(model::HydrogenAtom, r)`

```math
\begin{aligned}
  V(r)
  &= - \frac{Ze^2}{4\pi\varepsilon_0 r} 
  &= - \frac{e^2}{4\pi\varepsilon_0 a_0} \frac{Z}{r/a_0}
  &= - \frac{Z}{r/a_0} E_\mathrm{h},
\end{aligned}
```
The domain is $0\leq r \lt \infty$.
""" V(model::HydrogenAtom, r)

@doc raw"""
`E(model::HydrogenAtom; n=1)`

```math
E_n
= -\frac{m_\mathrm{e} e^4 Z^2}{2n^2(4\pi\varepsilon_0)^2\hbar^2}
= -\frac{Z^2}{2n^2} E_\mathrm{h},
```
where ``E_\mathrm{h}`` is the Hartree energy, one of atomic unit. About atomic units, see section 3.9.2 of the [IUPAC GreenBook](https://iupac.org/what-we-do/books/greenbook/). In other units, ``E_\mathrm{h} = 27.211~386~245~988(53)~\mathrm{eV}`` from [here](https://physics.nist.gov/cgi-bin/cuu/Value?hrev).
""" E(model::HydrogenAtom; n=1)

@doc raw"""
`ψ(model::HydrogenAtom, r, θ, φ; n=1, l=0, m=0)`

```math
\psi_{nlm}(\pmb{r}) = R_{nl}(r) Y_{lm}(\theta,\varphi)
```
The domain is $0\leq r \lt \infty, 0\leq \theta \lt \pi, 0\leq \varphi \lt 2\pi$.
""" ψ(model::HydrogenAtom, r, θ, φ; n=1, l=0, m=0)

@doc raw"""
`R(model::HydrogenAtom, r; n=1, l=0)`

```math
R_{nl}(r) = -\sqrt{\frac{(n-l-1)!}{2n(n+l)!} \left(\frac{2Z}{n a_0}\right)^3} \left(\frac{2Zr}{n a_0}\right)^l \exp \left(-\frac{Zr}{n a_0}\right) L_{n+l}^{2l+1} \left(\frac{2Zr}{n a_0}\right),
```
where Laguerre polynomials are defined as ``L_n(x) = \frac{1}{n!} \mathrm{e}^x \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left( \mathrm{e}^{-x} x^n \right)``, and associated Laguerre polynomials are defined as ``L_n^{k}(x) = \frac{\mathrm{d}^k}{\mathrm{d}x^k} L_n(x)``. Note that replace ``2n(n+l)!`` with ``2n[(n+l)!]^3`` if Laguerre polynomials are defined as ``L_n(x) = \mathrm{e}^x \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left( \mathrm{e}^{-x} x^n \right)``. 
The domain is $0\leq r \lt \infty$.
""" R(model::HydrogenAtom, r; n=1, l=0)

@doc raw"""
`L(model::HydrogenAtom, x; n=0, k=0)`

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
""" L(model::HydrogenAtom, x; n=0, k=0)

@doc raw"""
`Y(model::HydrogenAtom, θ, φ; l=0, m=0)`

```math
Y_{lm}(\theta,\varphi) = (-1)^{\frac{|m|+m}{2}} \sqrt{\frac{2l+1}{4\pi} \frac{(l-|m|)!}{(l+|m|)!}} P_l^{|m|} (\cos\theta) \mathrm{e}^{im\varphi}.
```
The domain is $0\leq \theta \lt \pi, 0\leq \varphi \lt 2\pi$. Note that some variants are connected by 
```math
i^{|m|+m} \sqrt{\frac{(l-|m|)!}{(l+|m|)!}} P_l^{|m|} = (-1)^{\frac{|m|+m}{2}} \sqrt{\frac{(l-|m|)!}{(l+|m|)!}} P_l^{|m|} = (-1)^m \sqrt{\frac{(l-m)!}{(l+m)!}} P_l^{m}.
```

""" Y(model::HydrogenAtom, θ, φ; l=0, m=0)

@doc raw"""
`P(model::HydrogenAtom, x; n=0, m=0)`

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
""" P(model::HydrogenAtom, x; n=0, m=0)