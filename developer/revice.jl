# run `include("./developer/revice.jl")`
using Pkg
Pkg.instantiate()
# Pkg.add("Revise")
using Revise
Pkg.activate("./")
using Antique