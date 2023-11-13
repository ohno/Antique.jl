using Antiq
using Documenter

DocMeta.setdocmeta!(Antiq, :DocTestSetup, :(using Antiq); recursive=true)

makedocs(;
  modules=[Antiq],
  authors="Shuhei Ohno",
  repo="https://github.com/ohno/Antiq.jl/blob/{commit}{path}#{line}",
  sitename="Antiq.jl",
  format=Documenter.HTML(;
    prettyurls=get(ENV, "CI", "false") == "true",
    canonical="https://ohno.github.io/Antiq.jl",
    edit_link="main",
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
  ],
)

deploydocs(;
  repo="github.com/ohno/Antiq.jl",
  devbranch="main",
)
