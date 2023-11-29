module Antiq

  export antiq

  # Update this list when you add a model.
  models = [
    :InfinitePotentialWell,
    :HarmonicOscillator,
    :MorsePotential,
    :HydrogenAtom,
  ]

  # include statements
  for model in models
    include("./$(model).jl")
  end

  # I/O function
  function save(path, text)
    mkpath(dirname(path))
    file = open(path, "w")
    Base.write(file, text)
    close(file)
  end

  # I/O function
  function load(path)
    file = open(path, "r")
    text = Base.read(file, String)
    close(file)
    return text
  end

  # main functions
  function antiq(model; parameters...)
    # check existence of model
    if model ∉ models
      throw(ErrorException("\`:$(model)\` is not in the supported models $(models)."))
    end
    # load source code file
    if Sys.iswindows()
      source = load("$(dirname(@__FILE__()))\\\\$(model).jl")
    else
      source = load("$(dirname(@__FILE__()))/$(model).jl")
    end
    # replace name of module
    for m in eachmatch(Regex("module $(model)"), source)
      source = replace(source, m.match => "module $(model)$(hash(parameters))")
    end
    # replace parameters
    for parameter in parameters
      for m in eachmatch(Regex("$(parameter[1]) =.*# change here!"), source)
        source = replace(source, m.match => "$(parameter[1]) = $(parameter[2])")
      end
      # println("$(parameter[1]) = $(parameter[2])")
    end
    # modules
    return Meta.eval(Meta.parse(source))
  end

end

