module Antiq

  export InfinitePotentialWell, HarmonicOscillator, MorsePotential, HydrogenAtom, IPW, HO, MP, HA

  include("./InfinitePotentialWell.jl")
  include("./HarmonicOscillator.jl")
  include("./MorsePotential.jl")
  include("./HydrogenAtom.jl")

  IPW = InfinitePotentialWell
  HO = HarmonicOscillator
  MP = MorsePotential
  HA = HydrogenAtom

end

