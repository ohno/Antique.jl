```@meta
CurrentModule = Antique
```

# Delta Potential

The Delta potential is one of the simplest models for quantum mechanical system in 1D.
It always has one bound state and its wave function has a cusp at the origin.

## Definitions

This model is described with the time-independent Schrödinger equation
```math
  \hat{H} \psi(x) = E \psi(x),
```
and the Hamiltonian
```math
  \hat{H} = - \frac{\hbar^2}{2m} \frac{\mathrm{d}^2}{\mathrm{d}x ^2} + V(x).
```
Parameters are specified with the following struct.

#### Parameters
```@docs; canonical=false
Antique.DeltaPotential
```

#### Potential
```@docs; canonical=false
Antique.V(::DeltaPotential, ::Any)
```

#### Eigen Values
```@docs; canonical=false
Antique.E(::DeltaPotential)
```

#### Eigen Functions
```@docs; canonical=false
Antique.ψ(::DeltaPotential, ::Any)
```

## Usage & Examples

[Install Antique.jl](@ref Install) for the first use and run `using Antique` before each use. The energy `E()`, wavefunction `ψ()`, potential `V()` and some other functions are suppoted. In this system, the model is generated by `DeltaPotential` and several parameters `α`, `m` and `ℏ` are set as optional arguments.

```julia
using Antique
DP = DeltaPotential(α=1.0, m=1.0, ℏ=1.0)
```




Parameters:

```julia
julia> DP.α
1.0

julia> DP.m
1.0

julia> DP.ℏ
1.0
```



Eigen values:

```julia
julia> E(DP)
-0.5
```



Wave functions:

```julia
DP = DeltaPotential(α=0.1, m=0.5, ℏ=0.1)
x = LinRange(-2,2,500);

using Plots
plot(x, x->ψ(DP,x), linewidth=3)
plot!(xlim=[-2,2], ylim=[0,2.5], legend=false)
plot!(xlabel="x", ylabel="ψ(x)", title="Delta Potential")
```

![](./assets/fig//DeltaPotential_4_1.png)



Potential energy curve, Energy levels, Wave functions:

```julia
DP = DeltaPotential(α=0.1, m=0.5, ℏ=0.1)
x = LinRange(-2,2,500);

using Plots
plot(xlim=[-2,2], ylim=[-1,2.0], legend=false, xlabel="\$x\$", ylabel="\$V(x),~E,~\\psi(x)+E\$")
plot!([-2,0,0,0,2], [0,0,-1,0,0], lw=1, lc=:black) # plot!(x, x->V(DP,x), lw=1, lc=:black)
plot!(x, x->ψ(DP,x) + E(DP), lw=2, lc=1)
hline!([E(DP)], lw=1, ls=:dash, lc=:black)
```

![](./assets/fig//DeltaPotential_5_1.png)



## Testing

Unit testing and Integration testing were done using numerical integration ([QuadGK.jl](https://juliamath.github.io/QuadGK.jl/stable/)). The test script is [here](https://github.com/ohno/Antique.jl/blob/main/test/DeltaPotential.jl).

#### Normalization of $\psi(x)$

```math
\int_{-\infty}^{\infty} \psi^\ast(x) \psi(x) ~\mathrm{d}x = 1
```

```
  α |   m |   ℏ |        analytical |         numerical 
--- | --- | --- | ----------------- | ----------------- 
0.1 | 0.1 | 0.1 |    1.000000000000 |    1.000000000000 ✔
0.1 | 0.1 | 1.0 |    1.000000000000 |    1.000000000000 ✔
0.1 | 0.1 | 7.0 |    1.000000000000 |    1.000004676239 ✔
0.1 | 1.0 | 0.1 |    1.000000000000 |    1.000000000000 ✔
0.1 | 1.0 | 1.0 |    1.000000000000 |    1.000000000000 ✔
0.1 | 1.0 | 7.0 |    1.000000000000 |    0.999999999999 ✔
0.1 | 7.0 | 0.1 |    1.000000000000 |    1.000000000000 ✔
0.1 | 7.0 | 1.0 |    1.000000000000 |    1.000000000000 ✔
0.1 | 7.0 | 7.0 |    1.000000000000 |    1.000000000000 ✔
1.0 | 0.1 | 0.1 |    1.000000000000 |    1.000000000000 ✔
1.0 | 0.1 | 1.0 |    1.000000000000 |    1.000000000000 ✔
1.0 | 0.1 | 7.0 |    1.000000000000 |    0.999999999999 ✔
1.0 | 1.0 | 0.1 |    1.000000000000 |    1.000000000000 ✔
1.0 | 1.0 | 1.0 |    1.000000000000 |    1.000000000000 ✔
1.0 | 1.0 | 7.0 |    1.000000000000 |    1.000000000000 ✔
1.0 | 7.0 | 0.1 |    1.000000000000 |    1.000000000000 ✔
1.0 | 7.0 | 1.0 |    1.000000000000 |    1.000000000000 ✔
1.0 | 7.0 | 7.0 |    1.000000000000 |    1.000000000000 ✔
7.0 | 0.1 | 0.1 |    1.000000000000 |    1.000000000000 ✔
7.0 | 0.1 | 1.0 |    1.000000000000 |    1.000000000000 ✔
7.0 | 0.1 | 7.0 |    1.000000000000 |    1.000000000000 ✔
7.0 | 1.0 | 0.1 |    1.000000000000 |    1.000000000000 ✔
7.0 | 1.0 | 1.0 |    1.000000000000 |    1.000000000000 ✔
7.0 | 1.0 | 7.0 |    1.000000000000 |    1.000000000000 ✔
7.0 | 7.0 | 0.1 |    1.000000000000 |    1.000000000000 ✔
7.0 | 7.0 | 1.0 |    1.000000000000 |    1.000000000000 ✔
7.0 | 7.0 | 7.0 |    1.000000000000 |    1.000000000000 ✔

```
