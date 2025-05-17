# Please run `include("./dev/docs.jl")` on RELP.

# run(`julia --project=docs/ -e 'cd("docs"); using Pkg; Pkg.activate("./"); Pkg.instantiate()'`) # Please try this command after removing Manifest.toml if you have any error.
run(`julia --project=docs/ -e 'cd("docs"); using Pkg; Pkg.activate("./"); include("make.jl")'`)
