```@meta
CurrentModule = Antique
```

# Antique.jl

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

```@raw html
<div class="catalog">
  <div class="item">
    <a target="_blank" href="./InfinitePotentialWell">
      <img src="assets/fig/InfinitePotentialWell_6_1.png" alt="InfinitePotentialWell"/>
    </a>
    <code>:InfinitePotentialWell</code>
  </div>
  <div class="item">
    <a target="_blank" href="./HarmonicOscillator">
      <img src="assets/fig/HarmonicOscillator_6_1.png" alt="HarmonicOscillator"/>
    </a>
    <code>:HarmonicOscillator</code>
  </div>
  <div class="item">
    <a target="_blank" href="./MorsePotential">
      <img src="assets/fig/MorsePotential_6_1.png" alt="MorsePotential"/>
    </a>
    <code>:MorsePotential</code>
  </div>
  <div class="item">
    <a target="_blank" href="./HydrogenAtom">
      <img src="assets/fig/HydrogenAtom_5_1.png" alt="HydrogenAtom"/>
    </a>
    <code>:HydrogenAtom</code>
  </div>
</div>
```

## Future Works

[List of quantum-mechanical systems with analytical solutions](https://en.wikipedia.org/wiki/List_of_quantum-mechanical_systems_with_analytical_solutions)

## Acknowledgment

This package was named by [@KB-satou](https://github.com/KB-satou) and [@ultimatile](https://github.com/ultimatile).