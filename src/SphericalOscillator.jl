export SphericalOscillator, V, E, ψ

@kwdef struct SphericalOscillator
  k = 1.0
  μ = 1.0
  ℏ = 1.0
end

function V(model::SphericalOscillator, r)
  if !(0 ≤ r)
    throw(DomainError("r = $r", "r must be non-negative: 0 ≤ r."))
  end
  k = model.k
  return 1/2 * k * r^2
end

function E(model::SphericalOscillator; n::Int=0, l::Int=0)
  if !(0 ≤ n && 0 ≤ l)
    throw(DomainError("(n,l) = ($n,$l)", "This function is defined for 0 ≤ n and 0 ≤ l"))
  end
  ℏ = model.ℏ
  μ = model.μ
  k = model.k 
  ω = sqrt(k/μ)
  return (2*n + l + 3/2) * ℏ * ω
end

function ψ(model::SphericalOscillator, r, θ, φ; n::Int=0, l::Int=0, m::Int=0)
  if !(0 ≤ n && 0 ≤ l && -l ≤ m ≤ l)
    throw(DomainError("(n,l,m) = ($n,$l,$m)", "This function is defined for 0 ≤ n, 0 ≤ l and -l ≤ m ≤ l."))
  end
  if !(0 ≤ r && 0 ≤ θ < π && 0 ≤ φ < 2π)
    throw(DomainError("(r,θ,φ) = ($r,$θ,$φ)", "This function is defined for 0 ≤ r, 0 ≤ θ < π, 0 ≤ φ < 2π."))
  end
  return R(model, r, n=n, l=l) * Y(model, θ, φ, l=l, m=m)
end

function R(model::SphericalOscillator, r; n=0, l=0)
  ℏ = model.ℏ
  μ = model.μ
  k = model.k 
  ω = sqrt(k/μ)
  γ = μ*ω/ℏ
  ξ = sqrt(γ)*abs(r)
  fact(n) = n>0 ? n*fact(n-1) : 1    # n!
  ffact(n) = n>0 ? n*ffact(n-2) : 1  # n!!
  N = sqrt( γ^(3/2)/(2*sqrt(π))) * sqrt( 2^(n+l+3) * fact(n)/ffact(2n+2l+1))
  return N * ξ^l * exp(-ξ^2/2) * L(model, ξ^2, n=n, α=l+1/2)
end

function L(model::SphericalOscillator, x; n=0, α=0)
  if isinteger(α)
    return sum(k -> (-1)^(k) * (Int(gamma(α+n+1)) // Int((gamma(α+1+k)*gamma(n-k+1)))) * x^k * 1 // factorial(k), 0:n)
  else
    return sum(k -> (-1)^(k) * (gamma(α+n+1) / (gamma(α+1+k)*gamma(n-k+1))) * x^k / factorial(k), 0:n)
  end
end

function Y(model::SphericalOscillator, θ, φ; l=0, m=0)
  N = (im)^(m+abs(m)) * sqrt( (2*l+1)*factorial(l-Int(abs(m))) / (2*factorial(l+Int(abs(m)))) )
  return N * P(model,cos(θ), n=l, m=Int(abs(m))) * exp(im*m*φ) / sqrt(2*π)
end

function P(model::SphericalOscillator, x; n=0, m=0)
  return (1//2)^n * (1-x^2)^(m//2) * sum(j -> (-1)^j * factorial(2*n-2*j) // (factorial(j) * factorial(n-j) * factorial(n-2*j-m)) * x^(n-2*j-m), 0:Int(floor((n-m)/2)))
end


@doc raw"""
`SphericalOscillator(k=1.0, μ=1.0, ℏ=1.0)`

``k`` is the force constant, ``μ`` is the mass of particle and ``\hbar`` is the reduced Planck constant (Dirac's constant).
""" SphericalOscillator

@doc raw"""
`V(model::SphericalOscillator, r)`

```math
V(r)
= \frac{1}{2} k r^2
= \frac{1}{2} \mu \omega^2 r^2
= \frac{1}{2} \hbar \omega \xi^2,
```
where ``\omega = \sqrt{k/\mu}`` is the angular frequency and ``\xi = \sqrt{\frac{\mu\omega}{\hbar}}r``.
""" V(model::SphericalOscillator, x)

@doc raw"""
`E(model::SphericalOscillator; n=0, l=0)`

```math
E_{nl}
= \left(2n + l + \frac{3}{2}\right)\hbar \omega,
```
where ``\omega = \sqrt{k/\mu}``.
""" E(model::SphericalOscillator; n=0, l=0)

@doc raw"""
`ψ(model::SphericalOscillator, r, θ, φ; n=0, l=0, m=0)`

```math
\psi_{nlm}(\pmb{r}) = R_{nl}(r) Y_{lm}(\theta,\varphi)
```
The domain is $0\leq r \lt \infty, 0\leq \theta \lt \pi, 0\leq \varphi \lt 2\pi$.
""" ψ(model::SphericalOscillator, r, θ, φ; n=0, l=0, m=0)

@doc raw"""
`R(model::SphericalOscillator, r; n=0, l=0)`

```math
R_{nl}(r) = \sqrt{ \frac{\gamma^{3/2}}{2\sqrt{\pi}}} \sqrt{\frac{2^{n+l+3} n!}{(2n+2l+1)!!}} \xi^l \exp\left(-\xi^2/2\right)L_{n}^{(l+\frac{1}{2})} \left(\xi^2\right),
```
where ``\gamma = \mu\omega/\hbar`` and ``\xi = \sqrt{\gamma}r = \sqrt{\mu\omega/\hbar}r`` are defined. The generalized Laguerre polynomials are defined as ``L_n^{(\alpha)}(x) = \frac{x^{-\alpha} \mathrm{e}^x}{n !} \frac{\mathrm{d}^n}{\mathrm{d} x^n}\left(\mathrm{e}^{-x} x^{n+\alpha}\right)``. The domain is $0\leq r \lt \infty$.
""" R(model::SphericalOscillator, r; n=0, l=0)

@doc raw"""
`L(model::SphericalOscillator, x; n=0, α=0)`

!!! note
    The generalized Laguerre polynomials $L_n^{(\alpha)}(x)$, not the associated Laguerre polynomials $L_n^{k}(x)$, are used in this model.

Rodrigues' formula & closed-form:
```math
\begin{aligned}
  L_n^{(\alpha)}(x)
  &= \frac{x^{-\alpha}e^x}{n!} \frac{d^n}{dx^n}\left(x^{n+\alpha}e^{-x}\right) \\
  &= \sum_{k=0}^n(-1)^k \left(\begin{array}{l} n+\alpha \\ n-k \end{array}\right) \frac{x^k}{k !} \\
  &= \sum_{k=0}^n(-1)^k \frac{\Gamma(\alpha+n+1)}{\Gamma(\alpha+k+1)\Gamma(n-k+1)} \frac{x^k}{k !}.
\end{aligned}
```
Examples:
```math
\begin{aligned}
  L_0^{(0)}(x) &= 1, \\
  L_1^{(0)}(x) &= 1 - x, \\
  L_1^{(1)}(x) &= 2 - x, \\
  L_2^{(0)}(x) &= 1 - 2 x + 1/2 x^{2}, \\
  L_2^{(1)}(x) &= 3 - 3 x + 1/2 x^{2}, \\
  L_2^{(2)}(x) &= 6 - 4 x + 1/2 x^{2}, \\
  L_3^{(0)}(x) &= 1 - 3 x + 3/2 x^{2} - 1/6 x^{3}, \\
  L_3^{(1)}(x) &= 4 - 6 x + 2 x^{2} - 1/6 x^{3}, \\
  L_3^{(2)}(x) &= 10 - 10 x + 5/2 x^{2} - 1/6 x^{3}, \\
  L_3^{(3)}(x) &= 20 - 15 x + 3 x^{2} - 1/6 x^{3}, \\
  L_4^{(0)}(x) &= 1 - 4 x + 3 x^{2} - 2/3 x^{3} + 1/24 x^{4}, \\
  L_4^{(1)}(x) &= 5 - 10 x + 5 x^{2} - 5/6 x^{3} + 1/24 x^{4}, \\
  L_4^{(2)}(x) &= 15 - 20 x + 15/2 x^{2} - 1 x^{3} + 1/24 x^{4}, \\
  L_4^{(3)}(x) &= 35 - 35 x + 21/2 x^{2} - 7/6 x^{3} + 1/24 x^{4}, \\
  L_4^{(4)}(x) &= 70 - 56 x + 14 x^{2} - 4/3 x^{3} + 1/24 x^{4}, \\
  \vdots
\end{aligned}
```
""" L(model::SphericalOscillator, x; n=0, α=0)

@doc raw"""
`Y(model::SphericalOscillator, θ, φ; l=0, m=0)`

```math
Y_{lm}(\theta,\varphi) = (-1)^{\frac{|m|+m}{2}} \sqrt{\frac{2l+1}{4\pi} \frac{(l-|m|)!}{(l+|m|)!}} P_l^{|m|} (\cos\theta) \mathrm{e}^{im\varphi}.
```
The domain is $0\leq \theta \lt \pi, 0\leq \varphi \lt 2\pi$. Note that some variants are connected by 
```math
i^{|m|+m} \sqrt{\frac{(l-|m|)!}{(l+|m|)!}} P_l^{|m|} = (-1)^{\frac{|m|+m}{2}} \sqrt{\frac{(l-|m|)!}{(l+|m|)!}} P_l^{|m|} = (-1)^m \sqrt{\frac{(l-m)!}{(l+m)!}} P_l^{m}.
```

""" Y(model::SphericalOscillator, θ, φ; l=0, m=0)

@doc raw"""
`P(model::SphericalOscillator, x; n=0, m=0)`

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
""" P(model::SphericalOscillator, x; n=0, m=0)