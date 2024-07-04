```@meta
CurrentModule = Antique
```

# Antique.jl

Self-contained, Well-Tested, Well-Documented **An**aly**ti**cal Solutions of **Qu**antum Mechanical **E**quations.

## Install

Run the following code on the REPL or Jupyter Notebook to install this package.

```julia
]add Antique
```

Or specify the version like `]add Antique@0.11.0` to install a specific version. The version of this package can be found at `]status Antique`.

## Usage & Examples

[Install Antique.jl](@ref Install) for the first use and run `using Antique` before each use.

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
- [Infinite PotentialWell 3D](https://ohno.github.io/Antique.jl/stable/InfinitePotentialWell3D/) `InfinitePotentialWell3D`
- [Spherical Oscillator](https://ohno.github.io/Antique.jl/stable/SphericalOscillator/) `SphericalOscillator`
- [Hydrogen Atom](https://ohno.github.io/Antique.jl/stable/HydrogenAtom/) `HydrogenAtom`
- [Coulomb 2-Body System](https://ohno.github.io/Antique.jl/stable/HydrogenAtom/) `CoulombTwoBody`

## Demonstration

This is an example of a variational calculation for the hydrogen atom based on [Thijssen(2007)](https://doi.org/10.1017/CBO9781139171397). We check the accuracy of the numerical solution by comparison with the analytical solution. Comparing wavefunctions can be difficult, but Antique.jl makes it easy. You can extend it to excited states ($n>1$) as well as the ground state ($n=1$). Thus, Antique.jl is useful for testing numerical methods. We hope many numerical methods to be developed using Antique.jl.

```@example demonstration
# calculations based on Thijssen(2007) https://doi.org/10.1017/CBO9781139171397
using LinearAlgebra
α = [13.00773, 1.962079, 0.444529, 0.1219492] 
nₘₐₓ = length(α)
S = [(pi/(α[i]+α[j]))^(3/2) for i=1:nₘₐₓ, j=1:nₘₐₓ]
T = [3*pi^(3/2)*α[i]*α[j]/(α[i]+α[j])^(5/2) for i=1:nₘₐₓ, j=1:nₘₐₓ]
V = [-2*pi/(α[i]+α[j]) for i=1:nₘₐₓ, j=1:nₘₐₓ]
H = T + V
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
f = Figure(size=(420,300), fontsize=11.5)
ax = Axis(f[1,1], xlabel=L"$r$", ylabel=L"$\psi(r,0,0)$", limits=(0,4,0,0.6), ylabelsize=16.5, xlabelsize=16.5)
l1 = lines!(ax, 0:0.01:10, r -> sum(C[:,1] .* exp.(-α*r^2)))
l2 = lines!(ax, 0:0.01:10, r -> real(Antique.ψ(HA,r,0,0)), color=:black, linestyle=:dash, label="Antique.jl")
axislegend(ax, [l1,l2], ["Numerical, Thijssen(2007)","Analytical, Antique.jl"], position=:rt)
f
save("assets/fig/demonstration.png", f) # hide
; # hide
```

![](assets/fig/demonstration.png)

## Future Works

The candidate models are listed on the Wikipedia page of [List of quantum-mechanical systems with analytical solutions](https://en.wikipedia.org/wiki/List_of_quantum-mechanical_systems_with_analytical_solutions). Please submit your requests and suggestions as [issues on GitHub](https://github.com/ohno/Antique.jl/issues).

## Developer's Guide

[Here](https://github.com/ohno/Antique.jl?tab=readme-ov-file#developers-guide) is the guideline for adding new models.

## Acknowledgment

Thanks to all contributors. This package was named by [@KB-satou](https://github.com/KB-satou) and [@ultimatile](https://github.com/ultimatile). [@MartinMikkelsen](https://github.com/MartinMikkelsen) contributed to writing docstrings. Special thanks to [@hyrodium](https://github.com/hyrodium) for his help with managing the documentation and advice on coding style. [@lhapp27](https://github.com/lhapp27) implemented 2 models, and [@ajarifi](https://github.com/ajarifi) implemented 3 models.
