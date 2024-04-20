using Antique
using Documenter
using CairoMakie

DocMeta.setdocmeta!(Antique, :DocTestSetup, :(using Antique); recursive=true)

makedocs(;
  modules=[Antique],
  authors="Shuhei Ohno",
  repo="https://github.com/ohno/Antique.jl/blob/{commit}{path}#{line}",
  sitename="Antique.jl",
  format=Documenter.HTML(;
    prettyurls=get(ENV, "CI", "false") == "true",
    canonical="https://ohno.github.io/Antique.jl",
    edit_link="main",
    repolink="https://github.com/ohno/Antique.jl",
    assets=String[
      "./assets/fig/logo.ico",
      "./assets/css/catalog.css"
    ],
  ),
  pages=[
    "Home" => "index.md",
    "Delta Potential"              => "DeltaPotential.md"         ,
    "Infinite Potential Well"      => "InfinitePotentialWell.md"  ,
    "Harmonic Oscillator"          => "HarmonicOscillator.md"     ,
    "Morse Potential"              => "MorsePotential.md"         ,
    "Pöschl-Teller Potential"      => "PoschlTeller.md"           ,
    "Rigid Rotor"                  => "RigidRotor.md"             ,
    "Hydrogen Atom"                => "HydrogenAtom.md"           ,
    "Coulomb 2-Body System"        => "CoulombTwoBody.md"         ,
    # "Infinite Potential Well 3D" => "InfinitePotentialWell3D.md",
    # "Harmonic Oscillator 3D"     => "HarmonicOscillator3D.md"   ,
    # "API reference"              => "API.md"                    ,
  ],
)

deploydocs(;
  repo="github.com/ohno/Antique.jl",
  devbranch="main",
)
