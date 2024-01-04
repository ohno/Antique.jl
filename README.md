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

Install Antique.jl for the first use and run `using Antique` before each use. The function `antique(model, parameters...)` returns a module that has `E`, `ψ`, `V` and some other functions. Here are examples in hydrogen-like atom. The analytical notation of energy (eigen value of the Hamiltonian) is written as

```math
E_n = -\frac{Z^2}{2n^2} E_\mathrm{h}.
```

Hydrogen atom has symbol $\mathrm{H}$ and atomic number 1 ($Z=1$). Therefore the ground state ($n=1$) energy is $-\frac{1}{2} E_\mathrm{h}$.

```julia
using Antique
H = antique(:HydrogenAtom, Z=1)
H.E(n=1)
# output> -0.5
```

Helium cation has symbol $\mathrm{He}^+$ and atomic number 2 ($Z=2$). Therefore the ground state ($n=1$) energy is $-2 E_\mathrm{h}$.

```julia
using Antique
He⁺ = antique(:HydrogenAtom, Z=2)
He⁺.E(n=1)
# output> -2.0
```

There are more examples on each model page.

## Supported Models

- [Infinite Potential Well](https://ohno.github.io/Antique.jl/dev/InfinitePotentialWell/) `:InfinitePotentialWell`
- [Harmonic Oscillator](https://ohno.github.io/Antique.jl/dev/HarmonicOscillator/) `:HarmonicOscillator`
- [Morse Potential](https://ohno.github.io/Antique.jl/dev/MorsePotential/) `:MorsePotential`
- [Hydrogen Atom](https://ohno.github.io/Antique.jl/dev/HydrogenAtom/) `:HydrogenAtom`

## Future Works

[List of quantum-mechanical systems with analytical solutions](https://en.wikipedia.org/wiki/List_of_quantum-mechanical_systems_with_analytical_solutions)

## Acknowledgment

This package was named by [@KB-satou](https://github.com/KB-satou) and [@ultimatile](https://github.com/ultimatile).
