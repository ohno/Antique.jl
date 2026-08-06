module Antique

  # for Julia 1.1
  import Base:@kwdef

  # Export public functions
  export energy
  export potential
  export wavefunction
  export radial_function
  export spherical_harmonic
  export hermite_polynomial
  export laguerre_polynomial
  export legendre_polynomial
  export n_max

  # Export all models
  export InfinitePotentialWell
  export HarmonicOscillator
  export MorsePotential
  export HydrogenAtom
  export DeltaPotential
  export PoschlTeller
  export SphericalOscillator
  export RigidRotor
  export InfinitePotentialWell3D
  export CoulombTwoBody

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
    hermite_polynomial(model, args...; kwargs...)

  Evaluate a Hermite polynomial for a given `model`.
  """
  function hermite_polynomial end

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

  Evaluate a Laguerre polynomial for a given `model`.
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

  # Define the display representation for all models.
  Base.show(io::IO, model::AbstractModel) = print(io, "Antique.", nameof(typeof(model)), "(", join(("$(name)=$(getfield(model, name))" for name in fieldnames(typeof(model))), ", "), ")")

end
