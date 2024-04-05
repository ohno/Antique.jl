# Antique.jl

[![Build Status](https://github.com/ohno/Antique.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ohno/Antique.jl/actions/workflows/CI.yml?query=branch%3Amain) 
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ohno.github.io/Antique.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ohno.github.io/Antique.jl/dev/)

Self-contained, Well-Tested, Well-Documented **An**aly**ti**cal Solutions of **Qu**antum Mechanical **E**quations.

## Install

To install this package, run the following code in your Jupyter Notebook:

```julia
using Pkg; Pkg.add("Antique")
```

## Usage & Examples

Install Antique.jl for the first use and run `using Antique` before each use. The energy `E()`, wavefunction `ψ()`, potential `V()` and some other functions are suppoted. Here are examples in hydrogen-like atom. The analytical notation of energy (eigen value of the Hamiltonian) is written as

```math
E_n = -\frac{Z^2}{2n^2} E_\mathrm{h}.
```

Hydrogen atom has symbol $\mathrm{H}$ and atomic number 1 ($Z=1$). Therefore the ground state ($n=1$) energy is $-\frac{1}{2} E_\mathrm{h}$.

```julia
using Antique
H = HydrogenAtom(Z=1)
E(H)
# output> -0.5
```

Helium cation has symbol $\mathrm{He}^+$ and atomic number 2 ($Z=2$). Therefore the ground state ($n=1$) energy is $-2 E_\mathrm{h}$.

```julia
using Antique
He⁺ = HydrogenAtom(Z=2)
E(He⁺)
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
- [Hydrogen Atom](https://ohno.github.io/Antique.jl/stable/HydrogenAtom/) `HydrogenAtom`

## Future Works

[List of quantum-mechanical systems with analytical solutions](https://en.wikipedia.org/wiki/List_of_quantum-mechanical_systems_with_analytical_solutions)

## Developer's Guide

This is the guideline for adding new models.

1. First, please submit a new issue [here](https://github.com/ohno/Antique.jl/issues). We need to find orthodox references (textbooks or papers, not Wikipedia) for the analytical solutions (eigenvalues and eigenfunctions) before the development. This will take more time than you think.
2. Fork [the repository](https://github.com/ohno/Antique.jl) on GitHub.
3. Clone the forked repository to your local machine by Git.
4. Please create 3 files:

| files | comments |
| --- | --- |
| `src/ModelName.jl` | Write the source codes and docstrings in this file. The most helpful examples are harmonic oscillators for one-dimensional systems and hydrogen atoms for three-dimensional systems. We recommend that you copy these files. First we need to create a structure `struct ModelName` with the same name as the model name (The best way is Find & Replace). Because the function names conflict, you must always give the struct `ModelName` as the fisrt argument to V, E, ψ and other functions. Multi-dispatch avoids conflicts. We recommend using Revice.jl while coding. Run `include("./developer/revice.jl")` on the REPL or use dev.ipynb.  |
| `test/ModelName.jl` | Write test code in this file. At a minimum, it is recommended to check the normalization and the orthogonality of wavefunction using QuadGK.jl. |
| `docs/src/ModelName.md` | Write documnetation in this file. Include at least the definition of the Hamiltonian and the analytical solutions (eigenvalues and eigenfunctions). Calls a docstring in the source code. |

5. Please rewrite 2 files:

| files | comments |
| - | - |
| `src/Antique.jl` | Add the new model name `:ModelName` to the `models = [...]` array in this file. `:` is required at the beginning. |
| `docs/make.jl` | Add the new model into `pages=[...]` in this file. |

6. Execute `include("./developer/test.jl")` to run tests. It will take few minutes to complete.
7. Execute `include("./developer/docs.jl")` to compile documents. HTML files (docs/build/*.html) will be generated. Please check them with Chrome or any other web browsers.
8. Commit and Push the codes.
9. Submit a pull request on GitHub.

## Acknowledgment

This package was named by [@KB-satou](https://github.com/KB-satou) and [@ultimatile](https://github.com/ultimatile).
