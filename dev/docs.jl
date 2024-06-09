# Please run `include("./dev/docs.jl")` on RELP.

# run(`julia --project=docs/ -e 'try; using CairoMakie; catch; using Pkg; Pkg.add("CairoMakie"); end;'`)
# run(`julia --project=docs/ -e 'try; using Documenter; catch; using Pkg; Pkg.add("Documenter"); end;'`)

# run(`julia --project=docs/ -e 'using Pkg; Pkg.activate("./"); Pkg.resolve();'`)
# run(`julia --project=docs/ -e 'using Pkg; Pkg.activate("./"); Pkg.instantiate();'`)
run(`julia --project=docs/ -e 'using Pkg; Pkg.activate("./"); cd("docs"); include("make.jl")'`)
