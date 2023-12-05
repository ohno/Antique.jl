# Antique.jl

[![Build Status](https://github.com/ohno/Antique.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ohno/Antique.jl/actions/workflows/CI.yml?query=branch%3Amain) 
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ohno.github.io/Antique.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ohno.github.io/Antique.jl/dev/)

Self-contained, Well-Tested, Well-Documented **An**aly**ti**cal Solutions of **Qu**antum Mechanical **E**quations.

## Install

To install this package, run the following code in your Jupyter Notebook:

```julia
using Pkg; Pkg.add(path="https://github.com/ohno/Antique.jl.git")
```

## Usage

To use this package, run the following code before each use:

```julia
using Antique
```

The function `antique(model, parameters...)` returns a module. Each module has `E()`, `ψ(x)` and some other functions.

## Examples

The energy of $1\mathrm{S}$ state in $\mathrm{H}$:
```julia
julia> H = antique(:HydrogenAtom, Z=1)
julia> H.E(n=1)
-0.5
```

The energy of $1\mathrm{S}$ state in $\mathrm{He}^+$:
```julia
julia> He⁺ = antique(:HydrogenAtom, Z=2)
julia> He⁺.E(n=1)
-2.0
```

## Supported Models

- [Infinite Potential Well](https://ohno.github.io/Antique.jl/dev/InfinitePotentialWell/) `:InfinitePotentialWell`
- [Harmonic Oscillator](https://ohno.github.io/Antique.jl/dev/HarmonicOscillator/) `:HarmonicOscillator`
- [Morse Potential](https://ohno.github.io/Antique.jl/dev/MorsePotential/) `:MorsePotential`
- [Hydrogen Atom](https://ohno.github.io/Antique.jl/dev/HydrogenAtom/) `:HydrogenAtom`

## Future Works

[List of quantum-mechanical systems with analytical solutions](https://en.wikipedia.org/wiki/List_of_quantum-mechanical_systems_with_analytical_solutions)

## Acknowledgment

This package was named by [@KB-satou](https://github.com/KB-satou) and [@ultimatile](https://github.com/ultimatile).
