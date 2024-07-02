export RigidRotor, V, E, ψ

# parameters
@kwdef struct RigidRotor
  m₁ = 1.0
  m₂ = 1.0
  R  = 1.0
  ℏ  = 1.0
end

# potential
function V(model::RigidRotor, r)
  if !(0 ≤ r)
    throw(DomainError("r = $r", "r must be non-negative: 0 ≤ r."))
  end
  return 0
end

# eigenvalues
function E(model::RigidRotor; l::Int=0)
  if !(0 ≤ l)
    throw(DomainError("l = $l", "l must be non-negative: 0 ≤ l."))
  end
  m₁ = model.m₁
  m₂ = model.m₂
  R = model.R
  ℏ = model.ℏ
  μ = (1/m₁ + 1/m₂)^(-1)
  I = μ * R^2
  return ℏ^2/(2*I) *l*(l+1)
end

# eigenfunctions
function ψ(model::RigidRotor, θ, φ; l::Int=0, m::Int=0)
  if !(0 ≤ l && -l ≤ m ≤ l)
    throw(DomainError("(l,m) = ($l,$m)", "This function is defined for 0 ≤ l and -l ≤ m ≤ l."))
  end
  if !(0 ≤ θ < π && 0 ≤ φ < 2π)
    throw(DomainError("(θ,φ) = ($θ,$φ)", "This function is defined for 0 ≤ θ < π and 0 ≤ φ < 2π."))
  end
  return Y(model, θ, φ; l=l, m=m)
end

# spherical harmonics
function Y(model::RigidRotor, θ, φ; l=0, m=0)
  N = (-1)^((abs(m)+m)/2) * sqrt( (2*l+1)*factorial(l-Int(abs(m))) / (2*factorial(l+Int(abs(m)))) )
  return N * P(model,cos(θ), n=l, m=Int(abs(m))) * exp(im*m*φ) / sqrt(2*π)
end

# associated Legendre polynomials
function P(model::RigidRotor, x; n=0, m=0)
  return (1//2)^n * (1-x^2)^(m//2) * sum(j -> (-1)^j * factorial(2*n-2*j) // (factorial(j) * factorial(n-j) * factorial(n-2*j-m)) * x^(n-2*j-m), 0:Int(floor((n-m)/2)))
end

# docstrings

@doc raw"""
This model is described with the time-independent Schrödinger equation
```math
  \hat{H} \psi(\theta,\varphi) = E \psi(\theta,\varphi),
```
and the Hamiltonian
```math
\begin{aligned}
  \hat{H} &= - \frac{\hbar^2}{2\mu} \nabla^2 + V(r), \\
          &= - \frac{\hbar^2}{2I} \left[ \frac{1}{\sin\theta} \frac{\partial}{\partial\theta} \left(\sin\theta \frac{\partial}{\partial\theta}\right) + \frac{1}{\sin^2\theta} \frac{\partial^2}{\partial\phi^2}  \right], \\
          &= \frac{L^2}{2I}
\end{aligned}
```
where $I=\mu R^2$ is the moment of intertia, $\mu=\left(\frac{1}{m_1}+\frac{1}{m_2}\right)^{-1}$ is the reduced mass of two particles, $R$ is the distance between the two particles, and $L^2$ is the angular momentum operator. Parameters are specified with the following struct:

```
RR = RigidRotor(m₁=1.0, m₂=1.0, R=1.0, ℏ=1.0)
```

``m₁`` and ``m₂`` are mass of two particles, ``R`` is the distance, and ``\hbar`` is the reduced Planck constant (Dirac's constant).

References:
- [D. A. McQuarrie, J. D. Simon, _Physical chemistry : a molecular approach_ (University Science Books, 1997)](https://uscibooks.aip.org/books/physical-chemistry-a-molecular-approach/) p.173, 5.8 The Energy Levels of a Rigid Rotator Are ``E = \hbar^2 J(J+1) / 2I``
- [Chemistry Libre Texts](https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Physical_Chemistry_(LibreTexts)/05%3A_The_Harmonic_Oscillator_and_the_Rigid_Rotor/5.08%3A_The_Energy_Levels_of_a_Rigid_Rotor)
""" RigidRotor

@doc raw"""
`V(model::RigidRotor, r)`

```math
\begin{aligned}
  V(r)
  &= 0,
\end{aligned}
```
""" V(model::RigidRotor, r)

@doc raw"""
`E(model::RigidRotor; l::Int=0)`

```math
E_l
= \frac{\hbar^2}{2I}l(l+1),
```
where ``I=\mu R^2`` is the moment of inertia, ``R`` is the distance, and ``\mu`` is the reduced mass of the two particles.
""" E(model::RigidRotor; l::Int=0)

@doc raw"""
`ψ(model::RigidRotor, θ, φ; l::Int=0, m::Int=0)`

```math
\psi_{lm}(\theta,\varphi) = Y_{lm}(\theta,\varphi)
```
The wave function is the spherical harmonics. The domain is $0\leq \theta \lt \pi$ and $0\leq \varphi \lt 2\pi$.
""" ψ(model::RigidRotor, θ, φ; l::Int=0, m::Int=0)

@doc raw"""
`Y(model::RigidRotor, θ, φ; l=0, m=0)`

```math
Y_{lm}(\theta,\varphi) = (-1)^{\frac{|m|+m}{2}} \sqrt{\frac{2l+1}{4\pi} \frac{(l-|m|)!}{(l+|m|)!}} P_l^{|m|} (\cos\theta) \mathrm{e}^{im\varphi}.
```
The domain is $0\leq \theta \lt \pi, 0\leq \varphi \lt 2\pi$. Note that some variants are connected by 
```math
i^{|m|+m} \sqrt{\frac{(l-|m|)!}{(l+|m|)!}} P_l^{|m|} = (-1)^{\frac{|m|+m}{2}} \sqrt{\frac{(l-|m|)!}{(l+|m|)!}} P_l^{|m|} = (-1)^m \sqrt{\frac{(l-m)!}{(l+m)!}} P_l^{m}.
```

""" Y(model::RigidRotor, θ, φ; l=0, m=0)

@doc raw"""
`P(model::RigidRotor, x; n=0, m=0)`

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
""" P(model::RigidRotor, x; n=0, m=0)