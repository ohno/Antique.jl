# run `include("./developer/docs.jl")`
run(`julia --project=docs/ -e 'try; using Plots; catch; using Pkg; Pkg.add("Plots"); end;'`)
run(`julia --project=docs/ -e 'try; using Documenter; catch; using Pkg; Pkg.add("Documenter"); end;'`)
run(`julia --project=docs/ -e 'using Pkg; Pkg.activate("./"); cd("docs"); include("make.jl")'`)
