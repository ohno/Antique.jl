using Antique
using Documenter
using Plots

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
      "./assets/logo.ico",
      "./assets/css/catalog.css"
    ],
  ),
  pages=[
    "Home" => "index.md",
    "Infinite Potential Well" => "InfinitePotentialWell.md",
    "Harmonic Oscillator"     => "HarmonicOscillator.md"   ,
    "Morse Potential"         => "MorsePotential.md"       ,
    "Hydrogen Atom"           => "HydrogenAtom.md"         ,
    "Delta Potential"         => "DeltaPotential.md"       ,
    # "Pöschl-Teller Potential" => "PoschlTeller.md"         ,
    # "API reference"           => "API.md"                  ,
  ],
)

deploydocs(;
  repo="github.com/ohno/Antique.jl",
  devbranch="main",
)
