```@meta
CurrentModule = Antiq
```

# Antiq

Self-contained, Well-Tested, Well-Documented Functions of Quantum Mechanical Models

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

Each module has `E()`, `ψ(x)` and other functions. An example of the energy of $1S$ state in hydrogen atom:

```julia
julia> HydrogenAtom.E(n=1)
-0.5
```

## Modules

```@index
```
