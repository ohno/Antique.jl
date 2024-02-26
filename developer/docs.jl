# run `include("./developer/docs.jl")`
run(`julia --project=docs/ -e 'cd("Antique.jl"); using Pkg; Pkg.activate("./"); cd("docs"); include("jmd2md.jl")'`)
run(`julia --project=docs/ -e 'cd("Antique.jl"); using Pkg; Pkg.activate("./"); cd("docs"); include("make.jl")'`)
