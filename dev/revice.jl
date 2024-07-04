# Please run `include("./dev/revice.jl")` on RELP.

# using Pkg
# Pkg.add("Revise")
using Revise

dir = dirname(@__FILE__) * "/../"
cd(dir)
@show pwd()
using Pkg
Pkg.activate(dir)
using Antique
