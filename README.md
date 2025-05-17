# Antique.jl

[![Build Status](https://github.com/ohno/Antique.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ohno/Antique.jl/actions/workflows/CI.yml?query=branch%3Amain) 
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ohno.github.io/Antique.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ohno.github.io/Antique.jl/dev/)

Antique.jl provides self-contained, well-tested, and well-documented implementations of **an**aly**ti**cal solutions to solvable **qu**antum m**e**chanical models. Analytical solutions are the most reliable benchmarks for software testing in the development of numerical methods. In addition to testing numerical methods, this package is useful for teaching quantum mechanics. We aim to support researchers, lecturers, students, and any person who is interested in quantum mechanics.

## Install

Run the following code on the REPL to install this package.

```julia
]add Antique@0.12.0
```

Or run `import Pkg; Pkg.add(; name="Antique", version="0.12.0")` to install on Jupyter Notebook. The version of this package can be found at `]status Antique` or `import Pkg; Pkg.status("Antique")`.

## Usage & Examples

Install Antique.jl for the first use and run `using Antique` before each use.

```julia
using Antique
```

The energy `E()`, the wave function `ψ()`, the potential `V()` and some other functions will be exported. There are two ways to avoid function name conflicts. Run `import Antique` instead of `using Antique`, and use the energy `Antique.E()`, the wave function `Antique.ψ()` and the potential `Antique.V()`. Or try giving other function names like `using Antique: V as potential, E as energy, ψ as wavefuntion, HydrogenAtom`. Here are examples for the hydrogen-like atom. The analytical notation of the energy (the eigen value of the Hamiltonian) is written as

```math
E_n = -\frac{Z^2}{2n^2} E_\mathrm{h}.
```

The Hydrogen atom has the symbol $\mathrm{H}$ and atomic number 1 ($Z=1$). Therefore the ground state ($n=1$) energy is $-\frac{1}{2} E_\mathrm{h}$.

```julia
H = HydrogenAtom(Z=1)
E(H, n=1)
# output> -0.5
```

The Helium cation has the symbol $\mathrm{He}^+$ and atomic number 2 ($Z=2$). Therefore the ground state ($n=1$) energy is $-2 E_\mathrm{h}$.

```julia
He⁺ = HydrogenAtom(Z=2)
E(He⁺, n=1)
# output> -2.0
```

There are more examples on each model page.

## Supported Models

- [Delta Potential](https://ohno.github.io/Antique.jl/stable/DeltaPotential/) `DeltaPotential`
- [Infinite Potential Well](https://ohno.github.io/Antique.jl/stable/InfinitePotentialWell/) `InfinitePotentialWell`
- [Harmonic Oscillator](https://ohno.github.io/Antique.jl/stable/HarmonicOscillator/) `HarmonicOscillator`
- [PoschlTeller](https://ohno.github.io/Antique.jl/stable/PoschlTeller/) `PoschlTeller`
- [Morse Potential](https://ohno.github.io/Antique.jl/stable/MorsePotential/) `MorsePotential`
- [Rigid Rotor](https://ohno.github.io/Antique.jl/stable/RigidRotor/) `RigidRotor`
- [Infinite PotentialWell 3D](https://ohno.github.io/Antique.jl/stable/InfinitePotentialWell3D/) `InfinitePotentialWell3D`
- [Spherical Oscillator](https://ohno.github.io/Antique.jl/stable/SphericalOscillator/) `SphericalOscillator`
- [Hydrogen Atom](https://ohno.github.io/Antique.jl/stable/HydrogenAtom/) `HydrogenAtom`
- [Coulomb 2-Body System](https://ohno.github.io/Antique.jl/stable/HydrogenAtom/) `CoulombTwoBody`

## Demonstration

This is an example of a variational calculation for the hydrogen atom based on [Thijssen(2007)](https://doi.org/10.1017/CBO9781139171397). We check the accuracy of the numerical solution by comparison with the analytical solution. Comparing wavefunctions can be difficult, but Antique.jl makes it easy. You can extend it to excited states ($n>1$) as well as the ground state ($n=1$). Thus, Antique.jl is useful for testing numerical methods. We hope many numerical methods to be developed using Antique.jl.

```julia
# calculations based on Thijssen(2007) https://doi.org/10.1017/CBO9781139171397
using LinearAlgebra
α = [13.00773, 1.962079, 0.444529, 0.1219492] 
nₘₐₓ = length(α)
S = [(pi/(α[i]+α[j]))^(3/2) for i=1:nₘₐₓ, j=1:nₘₐₓ]
H = [3*pi^(3/2)*α[i]*α[j]/(α[i]+α[j])^(5/2) - 2*pi/(α[i]+α[j]) for i=1:nₘₐₓ, j=1:nₘₐₓ]
E, C = eigen(Symmetric(H),Symmetric(S))

# norm & energy
import Antique
HA = Antique.HydrogenAtom(Z=1, Eₕ=1.0, a₀=1.0, mₑ=1.0, ℏ=1.0)
println("Norm")
println("  numerical : ", transpose(C[:,1]) * S * C[:,1])
println("  analytical: ", 1)
println("Energy")
println("  numerical : ", E[1])
println("  analytical: ", Antique.E(HA,n=1))

# wave function
using CairoMakie
fig = Figure(size=(420,300), fontsize=11, backgroundcolor=:transparent)
axis = Axis(fig[1,1], xlabel=L"$r$", ylabel=L"$\psi(r,0,0)$", limits=(0,4,0,0.6), ylabelsize=16.5, xlabelsize=16.5)
lines!(axis, 0:0.01:10, r -> sum(C[:,1] .* exp.(-α*r^2)), label="Numerical, Thijssen(2007)")
lines!(axis, 0:0.01:10, r -> real(Antique.ψ(HA,r,0,0)), color=:black, linestyle=:dash, label="Analytical, Antique.jl")
axislegend(axis, position=:rt, framevisible=false)
fig
```

```
Norm
  numerical : 0.9999999999999997
  analytical: 1
Energy
  numerical : -0.49927840566748566
  analytical: -0.5
```

![](docs/src/assets/fig/demonstration.png)

## Future Works

The candidate models are listed on the Wikipedia page of [List of quantum-mechanical systems with analytical solutions](https://en.wikipedia.org/wiki/List_of_quantum-mechanical_systems_with_analytical_solutions). Please submit your requests and suggestions as [issues on GitHub](https://github.com/ohno/Antique.jl/issues).

## Developer's Guide

This is the guideline for adding new models. Adding a new model may take from a few days to a few weeks due to reference search, test implementation, and writing documentation.

1. First, please submit a new issue or comment [here](https://github.com/ohno/Antique.jl/issues). I will assign you to the issue. We need to find orthodox references (textbooks or papers, not Wikipedia) for the analytical solutions (eigenvalues and eigenfunctions) before the development. This will take more time than you think.
2. Fork [the repository](https://github.com/ohno/Antique.jl) on GitHub.
3. Clone the forked repository to your local machine by Git.
4. Please create 3 files:

| files | comments |
| --- | --- |
| `src/ModelName.jl` | Write the source codes and docstrings in this file. The most helpful examples are the harmonic oscillator for one-dimensional systems and the hydrogen atom for three-dimensional systems. We recommend that you copy these files. First create a structure `struct ModelName` with the same name as the model name (The best way is Find & Replace). Because the function names conflict, you must always give the struct `ModelName` as the first argument to V, E, ψ and other functions. Multi-dispatch avoids conflicts. We recommend using Revice.jl while coding. Run `include("./dev/revice.jl")` on the REPL or use dev.ipynb. |
| `test/ModelName.jl` | Write test code in this file. At a minimum, please check the normalization and the orthogonality of the eigenfunctions using QuadGK.jl. Please also do tests for the eigenvalues (for example, calculate the expectation values of the Hamiltonian (energy) using the eigenfunctions and check that these values match the eigenvalues). |
| `docs/src/ModelName.md` | Write documentation in this file. Include at least the definition of the Hamiltonian and the analytical solutions (eigenvalues and eigenfunctions). Call a docstring in the source code (`src/ModelName.jl`) . |

5. Please rewrite 5 files:

| files | comments |
| - | - |
| `src/Antique.jl` | Add the new model name `:ModelName` to the `models = [...]` array in this file. `:` is required at the beginning. |
| `docs/make.jl` | Add the new model into `pages=[...]` in this file. |
| `test/runtests.jl` | Change `for model in [...]` in this file. Please test all models before pull requests. |
| `README.md` | Add the new model to the list of supported models. |
| `docs/index.md` | Add the new model to the list of supported models. |

6. Execute `include("./dev/test.jl")` to run tests. It will take few minutes to complete.
7. Execute `include("./dev/docs.jl")` to compile documents. HTML files (docs/build/*.html) will be generated. Please check them with Chrome or any other web browsers.
8. Commit and Push the codes.
9. Submit a pull request on GitHub.

## Acknowledgment

Thanks to all contributors. This package was named by [@KB-satou](https://github.com/KB-satou) and [@ultimatile](https://github.com/ultimatile). [@MartinMikkelsen](https://github.com/MartinMikkelsen) contributed to writing docstrings. Special thanks to [@hyrodium](https://github.com/hyrodium) for his help with managing the documentation and advice on coding style. [@lhapp27](https://github.com/lhapp27) implemented 2 models, and [@ajarifi](https://github.com/ajarifi) implemented 3 models.
