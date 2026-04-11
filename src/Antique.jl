module Antique

  # Update this list when you add a model.
  models = [
    :InfinitePotentialWell,
    :HarmonicOscillator,
    :MorsePotential,
    :HydrogenAtom,
    :DeltaPotential,
    :PoschlTeller,
    :SphericalOscillator,
    :RigidRotor,
    :InfinitePotentialWell3D,
    :CoulombTwoBody,
  ]

  # for Julia 1.1
  import Base:@kwdef

  # include statements
  for model in models
    include("./$(model).jl")
  end

  # override Base.string and Base.show
  for model in Antique.models
    Base.string(t::eval(Symbol(model))) = "Antique.$(typeof(t))(" * join(["$(symbol)=$(getproperty(t,symbol))" for symbol in fieldnames(typeof(t))], ", ") * ")"
    Base.show(io::IO, t::eval(Symbol(model))) = print(io, Base.string(t))
  end

end
