export HarmonicOscillator3D, V, E, ψ, R, L, Y, P

@kwdef struct HarmonicOscillator3D
    k = 1.0
    m = 1.0
    ℏ = 1.0
  end

function V(model::HarmonicOscillator3D, r)
  # if r<0
  #   throw(DomainError(r, "r=$r is out of the domain (0≦r)"))
  # end
  k = model.k
  return 1/2 * k * r^2
end

function E(model::HarmonicOscillator3D; n=1)
  ℏ = model.ℏ
  m = model.m
  k = model.k 
  ω = sqrt(k/m)
  return (n + 3/2) * ℏ * ω
end

function ψ(model::HarmonicOscillator3D, r, θ, φ; n=1, l=0, m=0)
  # if r<0
  #   throw(DomainError(r, "r=$r is out of the domain (0≦r)"))
  # end
  ℏ = model.ℏ
  μ = model.m #conflict with magnetic m
  k = model.k 
  return R(model, r, n=n, l=l) * Y(model, θ, φ, l=l, m=m)
end


function R(model::HarmonicOscillator3D, r; n=1, l=0)
  # if r<0
  #   throw(DomainError(r, "r=$r is out of the domain (0≦r)"))
  # end
  ℏ = model.ℏ
  m = model.m
  k = model.k 
  ω = sqrt(k/m)
  γ = m*ω/ℏ
  ξ = sqrt(γ)*abs(r)
  fact(n) = n>0 ? n*fact(n-1) : 1    # n!
  ffact(n) = n>0 ? n*ffact(n-2) : 1  # n!!
  N = sqrt( γ^(3/2)/(2*sqrt(π))) * sqrt( 2^(n+l+3) * fact(n)/ffact(2n+2l+1))
  return N * ξ^l * exp(-ξ^2/2) * L(model, ξ^2, n=n, k=l+1/2)
end

function L(model::HarmonicOscillator3D, x; n=0, k=0)
  if n < 0
      error("n must be a non-negative integer")
  end
  
  if n == 0
      return 1
  elseif n == 1
      return 1 + k - x
  else
      Lnm1 = 1
      Ln = 1 + k - x
      for m in 1:n-1
          Lnp1 = ((2*m + 1 + k - x) * Ln - (m + k) * Lnm1) / (m + 1)
          Lnm1, Ln = Ln, Lnp1
      end
      return Ln
  end
end 

function Y(model::HarmonicOscillator3D, θ, φ; l=0, m=0)
  N = (im)^(m+abs(m)) * sqrt( (2*l+1)*factorial(l-Int(abs(m))) / (2*factorial(l+Int(abs(m)))) )
  return N * P(model,cos(θ), n=l, m=Int(abs(m))) * exp(im*m*φ) / sqrt(2*π)
end

function P(model::HarmonicOscillator3D, x; n=0, m=0)
  return (1//2)^n * (1-x^2)^(m//2) * sum(j -> (-1)^j * factorial(2*n-2*j) // (factorial(j) * factorial(n-j) * factorial(n-2*j-m)) * x^(n-2*j-m), 0:Int(floor((n-m)/2)))
end


@doc raw"""
`HarmonicOscillator(k=1.0, m=1.0, ℏ=1.0)`

``k`` is the force constant, ``m`` is the mass of particle and ``\hbar`` is the reduced Planck constant (Dirac's constant).
""" HarmonicOscillator3D

@doc raw"""
`V(model::HarmonicOscillator3D, r)`

```math
V(r)
= \frac{1}{2} k r^2
= \frac{1}{2} m \omega^2 r^2
= \frac{1}{2} \hbar \omega \xi^2,
```
where ``\omega = \sqrt{k/m}`` is the angular frequency and ``\xi = \sqrt{\frac{m\omega}{\hbar}}r``.
""" V(model::HarmonicOscillator, x)

@doc raw"""
`E(model::HarmonicOscillator3D; n=1)`

```math
E_n
= \left(n + \frac{3}{2}\right)\hbar \omega,
```
where ``\omega = \sqrt{k/m}``
""" E(model::HarmonicOscillator3D; n=1)

@doc raw"""
`ψ(model::HarmonicOscillator3D, r, θ, φ; n=1, l=0, m=0)`

```math
\psi_{nlm}(\pmb{r}) = R_{nl}(r) Y_{lm}(\theta,\varphi)
```
The domain is $0\leq r \lt \infty, 0\leq \theta \lt \pi, 0\leq \varphi \lt 2\pi$.
""" ψ(model::HarmonicOscillator3D, r, θ, φ; n=1, l=0, m=0)

@doc raw"""
`R(model::HarmonicOscillator3D, r; n=1, l=0)`

```math
R_{nl}(r) = \sqrt{ \frac{\gamma^{3/2}}{(2\sqrt{\pi}}} \sqrt{\frac{2^{n+l+3} n!}{(2n+2l+1)!!}} \xi^l \exp\left(-\xi^2/2\right)L_{n}^{l+\frac{1}{2}} \left(\xi^2\right),
```
where Laguerre polynomials are defined as ``L_n(x) = \frac{1}{n!} \mathrm{e}^x \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left( \mathrm{e}^{-x} x^n \right)``, and associated Laguerre polynomials are defined as ``L_n^{k}(x) = \frac{\mathrm{d}^k}{\mathrm{d}x^k} L_n(x)``. Note that replace ``2n(n+l)!`` with ``2n[(n+l)!]^3`` if Laguerre polynomials are defined as ``L_n(x) = \mathrm{e}^x \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left( \mathrm{e}^{-x} x^n \right)``. 
The domain is $0\leq r \lt \infty$.
""" R(model::HarmonicOscillator3D, r; n=1, l=0)

@doc raw"""
`L(model::HarmonicOscillator3D, x; n=0, k=0)`

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
""" L(model::HarmonicOscillator3D, x; n=0, k=0)

@doc raw"""
`Y(model::HarmonicOscillator3D, θ, φ; l=0, m=0)`

```math
Y_{lm}(\theta,\varphi) = (-1)^{\frac{|m|+m}{2}} \sqrt{\frac{2l+1}{4\pi} \frac{(l-|m|)!}{(l+|m|)!}} P_l^{|m|} (\cos\theta) \mathrm{e}^{im\varphi}.
```
The domain is $0\leq \theta \lt \pi, 0\leq \varphi \lt 2\pi$. Note that some variants are connected by 
```math
i^{|m|+m} \sqrt{\frac{(l-|m|)!}{(l+|m|)!}} P_l^{|m|} = (-1)^{\frac{|m|+m}{2}} \sqrt{\frac{(l-|m|)!}{(l+|m|)!}} P_l^{|m|} = (-1)^m \sqrt{\frac{(l-m)!}{(l+m)!}} P_l^{m}.
```

""" Y(model::HarmonicOscillator3D, θ, φ; l=0, m=0)

@doc raw"""
`P(model::HarmonicOscillator3D, x; n=0, m=0)`

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
""" P(model::HarmonicOscillator3D, x; n=0, m=0)