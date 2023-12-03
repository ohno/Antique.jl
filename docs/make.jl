using AnalyticalSolutions
using Documenter

DocMeta.setdocmeta!(AnalyticalSolutions, :DocTestSetup, :(using AnalyticalSolutions); recursive=true)

makedocs(;
  modules=[AnalyticalSolutions],
  authors="Shuhei Ohno",
  repo="https://github.com/ohno/AnalyticalSolutions.jl/blob/{commit}{path}#{line}",
  sitename="AnalyticalSolutions.jl",
  format=Documenter.HTML(;
    prettyurls=get(ENV, "CI", "false") == "true",
    canonical="https://ohno.github.io/AnalyticalSolutions.jl",
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
  repo="github.com/ohno/AnalyticalSolutions.jl",
  devbranch="main",
)
