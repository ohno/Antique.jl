```@meta
CurrentModule = Antique
```

# Rigid rotor

The linear rigid rotor model can be used in quantum mechanics to predict the rotational energy of a diatomic molecule. The rotational energy depends on the moment of inertia for the system, $I$.

## Definitions

This model is described with the time-independent Schrödinger equation
```math
  \hat{H} \psi(\pmb{r}) = E \psi(\pmb{r}),
```
and the Hamiltonian
```math
\begin{aligned}
  \hat{H} &= - \frac{\hbar^2}{2\mu} \nabla^2 + V(r), \\
          &= - \frac{\hbar^2}{2I} \left[ \frac{1}{\sin\theta} \frac{\partial}{\partial\theta} \left(\sin\theta \frac{\partial}{\partial\theta}\right) + \frac{1}{\sin^2\theta} \frac{\partial^2}{\partial\phi^2}  \right]
\end{aligned}
```
where $I=\mu R^2$ is the moment of intertia, $\mu=\left(\frac{1}{m_1}+\frac{1}{m_2}\right)^{-1}$ is the reduced mass of two particles and $R$ is the distance between the two particles. Parameters are specified with the following struct.

#### Parameters
```@docs; canonical=false
Antique.RigidRotor
```

#### Potential
```@docs; canonical=false
Antique.V(::RigidRotor, ::Any)
```

#### Eigen Values
```@docs; canonical=false
Antique.E(::RigidRotor)
```

#### Eigen Functions
```@docs; canonical=false
Antique.ψ(::RigidRotor, ::Any, ::Any)
```

#### Spherical Harmonics
```@docs; canonical=false
Antique.Y(::RigidRotor, ::Any, ::Any)
```

#### Associated Legendre Polynomials
```@docs; canonical=false
Antique.P(::RigidRotor, ::Any)
```

## Usage & Examples

[Install Antique.jl](@ref Install) for the first use and run `using Antique` before each use. The energy `E()`, wavefunction `ψ()`, potential `V()` and some other functions are suppoted. In this system, the model is generated by `RigidRotor` and several parameters `m₁`, `m₂`, `R` and `ℏ` are set as optional arguments.

```@example RR
using Antique
RR = RigidRotor(m₁=1.0, m₂=1.0, R=1.0, ℏ=1.0)
; #hide
```

Parameters:

```@repl RR
RR.m₁
RR.m₂
RR.R
RR.ℏ
```

Eigen values:

```@repl RR
E(RR, l=0)
E(RR, l=1)
E(RR, l=2)
```

Wave functions:

```@repl RR
ψ(RR, 0, 0, l=2, m=1)
ψ(RR, π/4, 0, l=2, m=1)
ψ(RR, π/4, π/2, l=2, m=1)
```

```@example RR
using CairoMakie

f = Figure(size=(400,400))
ax = PolarAxis(f[1,1], title=L"$\theta\mapsto|\psi_{2,1}(\theta,0)|^2$", rticklabelsvisible=false)
lines!(ax, 0..2pi, θ->abs(ψ(RR,θ,0,l=2,m=1))^2, linewidth=2)

f
```

```@example RR
using CairoMakie

f = Figure(size=(400,400))
ax = PolarAxis(f[1,1], title=L"$\theta\mapsto|\psi_{lm}(\theta,0)|^2$", rticklabelsvisible=false)
l1 = lines!(ax, 0..2pi, θ->abs(ψ(RR,θ,0,l=0,m=0))^2, linewidth=2)
l2 = lines!(ax, 0..2pi, θ->abs(ψ(RR,θ,0,l=1,m=0))^2, linewidth=2)
l3 = lines!(ax, 0..2pi, θ->abs(ψ(RR,θ,0,l=1,m=1))^2, linewidth=2)
l4 = lines!(ax, 0..2pi, θ->abs(ψ(RR,θ,0,l=2,m=1))^2, linewidth=2)
Legend(f[2,1], [l1,l2,l3,l4], [L"(l,m)=(0,0)",L"(1,0)",L"(1,1)",L"(2,1)"], framevisible=false, orientation=:horizontal, tellwidth=false, tellheight=true)

f
save("assets/fig/RigidRotor.png", f) # hide
; # hide
```

![](assets/fig/RigidRotor.png)

## Testing

Unit testing and Integration testing were done using computer algebra system ([Symbolics.jl](https://symbolics.juliasymbolics.org/stable/)) and numerical integration ([QuadGK.jl](https://juliamath.github.io/QuadGK.jl/stable/)). The test script is [here](https://github.com/ohno/Antique.jl/blob/main/test/RigidRotor.jl).

```@eval
using Markdown
using Antique
Markdown.parse(Antique.load("../../test/result/RigidRotor.log"))
```