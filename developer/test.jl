# run `include("./developer/test.jl")`
using Pkg
Pkg.activate("./")
Pkg.test()