```@meta
CurrentModule = Antiq
```

# Antiq

**An**aly**ti**cal soulutions of Schrödinger e**q**uations, named by [@KB-satou](https://github.com/KB-satou) and [@ultimatile](https://github.com/ultimatile)

Self-contained, Well-Tested, Well-Documented Functions for Quantum Mechanical Models

## Install

To add this package, run the following code in your Jupyter Notebook:

```julia
using Pkg
Pkg.add(path="https://github.com/ohno/Antiq.jl.git")
```

or use `add` upon REPL(Package mode):

```julia
]
add https://github.com/ohno/Antiq.jl.git
```

## Usage

To use this package, run the following code in your Jupyter Notebook or code:

```julia
using Antiq
```

The function `antiq(model, parameters...)` returns a module. Each module has `E()`, `ψ(x)` and some other functions.

## Examples

The energy of $1S$ state in hydrogen atom:
```julia
julia> H = antiq(:HydrogenAtom, Z=1)
julia> H.E(n=1)
-0.5
```

The energy of $1S$ state in helium atom:
```
julia> He = antiq(:HydrogenAtom, Z=2)
julia> He.E(n=1)
-2.0
```

## Supported Models

```@raw html
<div class="catalog">
  <div class="item">
    <a target="_blank" href="./InfinitePotentialWell">
      <img src="assets/fig/InfinitePotentialWell_6_1.png" alt="InfinitePotentialWell"/>
    </a>
    <code>InfinitePotentialWell</code>
  </div>
  <div class="item">
    <a target="_blank" href="./HarmonicOscillator">
      <img src="assets/fig/HarmonicOscillator_6_1.png" alt="HarmonicOscillator"/>
    </a>
    <code>HarmonicOscillator</code>
  </div>
  <div class="item">
    <a target="_blank" href="./MorsePotential">
      <img src="assets/fig/MorsePotential_6_1.png" alt="MorsePotential"/>
    </a>
    <code>MorsePotential</code>
  </div>
  <div class="item">
    <a target="_blank" href="./HydrogenAtom">
      <img src="assets/fig/HydrogenAtom_5_1.png" alt="HydrogenAtom"/>
    </a>
    <code>HydrogenAtom</code>
  </div>
</div>
```

## Future Works

[List of quantum-mechanical systems with analytical solutions](https://en.wikipedia.org/wiki/List_of_quantum-mechanical_systems_with_analytical_solutions)
