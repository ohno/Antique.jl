module Antique

  export energy, potential, wavefunction, radial_function, laguerre_polynomial, spherical_harmonic, legendre_polynomial, laguerre_polynomial, n_max
  export InfinitePotentialWell, HarmonicOscillator, MorsePotential, HydrogenAtom, DeltaPotential
  export PoschlTeller, SphericalOscillator, RigidRotor, InfinitePotentialWell3D, CoulombTwoBody

  """
    energy(model; numbers..., parameters...)

  Return an eigenvalue for a given `model`.
  """
  function energy end

  """
    potential(model, args...)

  Return the potential for a given `model`.
  """
  function potential end

  """
    wavefunction(model, args...; kwargs...)

  Return a wavefunction for a given `model`.
  """
  function wavefunction end

  """
    radial_function(model, args...; kwargs...)

  Return a radial function for a given `model`.
  """
  function radial_function end

  """
    laguerre_polynomial(model, args...; kwargs...)

  Evaluate a Laguerre polynomial for a given `model`.
  """
  function laguerre_polynomial end

  """
    spherical_harmonic(model, args...; kwargs...)

  Evaluate a spherical harmonic for a given `model`.
  """
  function spherical_harmonic end

  """
    legendre_polynomial(model, args...; kwargs...)

  Evaluate a Rodrigues-form polynomial/function for a given `model`.
  """
  function legendre_polynomial end

  """
    laguerre_polynomial(model, args...; kwargs...)

  Evaluate a Hermite polynomial for a given `model`.
  """
  function laguerre_polynomial end

  """
    n_max(model)

  Return the maximum quantum number for a given `model`.
  """
  function n_max end

  abstract type AbstractModel end

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

  include("./InfinitePotentialWell.jl")
  include("./HarmonicOscillator.jl")
  include("./MorsePotential.jl")
  include("./HydrogenAtom.jl")
  include("./DeltaPotential.jl")
  include("./PoschlTeller.jl")
  include("./SphericalOscillator.jl")
  include("./RigidRotor.jl")
  include("./InfinitePotentialWell3D.jl")
  include("./CoulombTwoBody.jl")

  using .InfinitePotentialWells: InfinitePotentialWell
  using .HarmonicOscillators: HarmonicOscillator
  using .MorsePotentials: MorsePotential
  using .HydrogenAtoms: HydrogenAtom
  using .DeltaPotentials: DeltaPotential
  using .PoschlTellers: PoschlTeller
  using .SphericalOscillators: SphericalOscillator
  using .RigidRotors: RigidRotor
  using .InfinitePotentialWell3Ds: InfinitePotentialWell3D
  using .CoulombTwoBodies: CoulombTwoBody

  # override Base.string and Base.show
  for model in Antique.models
    Base.string(t::eval(Symbol(model))) = "Antique.$(typeof(t))(" * join(["$(symbol)=$(getproperty(t,symbol))" for symbol in fieldnames(typeof(t))], ", ") * ")"
    Base.show(io::IO, t::eval(Symbol(model))) = print(io, Base.string(t))
  end

end
