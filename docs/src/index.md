```@meta
CurrentModule = Antique
```

# Antique.jl

Self-contained, Well-Tested, Well-Documented **An**aly**ti**cal Solutions of **Qu**antum Mechanical **E**quations.

## Install

To install this package, run the following code in your Jupyter Notebook:

```julia
using Pkg; Pkg.add("Antique")
```

## Usage & Examples

[Install Antique.jl](@ref Install) for the first use and run `using Antique` before each use. The energy `E()`, wavefunction `ψ()`, potential `V()` and some other functions are suppoted. Here are examples in hydrogen-like atom. The analytical notation of energy (eigen value of the Hamiltonian) is written as

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
He⁺ = HydrogenAtom(Z=2)
E(He⁺)
# output> -2.0
```

There are more examples on each model page.

## Supported Models

```@raw html
<div class="catalog">
  <div class="item">
    <a target="_blank" href="./InfinitePotentialWell">
      <img src="assets/fig/InfinitePotentialWell.png" alt="InfinitePotentialWell"/>
    </a>
    <code>InfinitePotentialWell</code>
  </div>
  <div class="item">
    <a target="_blank" href="./HarmonicOscillator">
      <img src="assets/fig/HarmonicOscillator.png" alt="HarmonicOscillator"/>
    </a>
    <code>HarmonicOscillator</code>
  </div>
  <div class="item">
    <a target="_blank" href="./PoschlTeller">
      <img src="assets/fig/PoschlTeller.png" alt="PoschlTeller"/>
    </a>
    <code>PoschlTeller</code>
  </div>
  <div class="item">
    <a target="_blank" href="./MorsePotential">
      <img src="assets/fig/MorsePotential.png" alt="MorsePotential"/>
    </a>
    <code>MorsePotential</code>
  </div>
</div>
```

- [Delta Potential](https://ohno.github.io/Antique.jl/stable/DeltaPotential/) `DeltaPotential`
- [Infinite Potential Well](https://ohno.github.io/Antique.jl/stable/InfinitePotentialWell/) `InfinitePotentialWell`
- [Harmonic Oscillator](https://ohno.github.io/Antique.jl/stable/HarmonicOscillator/) `HarmonicOscillator`
- [PoschlTeller](https://ohno.github.io/Antique.jl/stable/PoschlTeller/) `PoschlTeller`
- [Morse Potential](https://ohno.github.io/Antique.jl/stable/MorsePotential/) `MorsePotential`
- [Rigid Rotor](https://ohno.github.io/Antique.jl/stable/RigidRotor/) `RigidRotor`
- [Spherical Oscillator](https://ohno.github.io/Antique.jl/stable/SphericalOscillator/) `SphericalOscillator`
- [Hydrogen Atom](https://ohno.github.io/Antique.jl/stable/HydrogenAtom/) `HydrogenAtom`
- [Coulomb 2-Body System](https://ohno.github.io/Antique.jl/stable/HydrogenAtom/) `CoulombTwoBody`

## Future Works

[List of quantum-mechanical systems with analytical solutions](https://en.wikipedia.org/wiki/List_of_quantum-mechanical_systems_with_analytical_solutions)

## Developer's Guide

[Here](https://github.com/ohno/Antique.jl?tab=readme-ov-file#developers-guide) is the guideline for adding new models.

## Acknowledgment

This package was named by [@KB-satou](https://github.com/KB-satou) and [@ultimatile](https://github.com/ultimatile).