```@meta
CurrentModule = Antique
```


# Delta Potential

The Delta potential is one the simplest models for quantum mechanical system in 1D.
It always has one bound state and its wave function has a cusp at the origin.

## Defitions

``\alpha`` is the potential strength, ``m`` is the mass of particle.

#### Schrödinger Equation
```math
  \hat{H} \psi(x) = E \psi(x)
```

#### Hamiltonian
```math
  \hat{H} = \frac{\hbar^2}{2m} \frac{\mathrm{d}^2}{\mathrm{d}x ^2} + V(x)
```

#### Potential
`V(model::DeltaPotential; x)`
```math
  V(x) = -\alpha \delta(x).
```

#### Eigen Values
`E(model::DeltaPotential)`
```math
  E = - \frac{m\alpha^2}{2\hbar^2}
```

#### Eigen Functions
`ψ(model::DeltaPotential, x)`
```math
   \psi(x) = \frac{\sqrt{m\alpha}}{\hbar} \mathrm{e}^{-m\alpha |x|/\hbar^2}
```

## Testing

Unit testing and Integration testing were done using numerical integration ([QuadGK.jl](https://juliamath.github.io/QuadGK.jl/stable/)). The test script is [here](https://github.com/ohno/Antique.jl/blob/main/test/DeltaPotential.jl).

Error: UndefVarError: `Antique` not defined
```


