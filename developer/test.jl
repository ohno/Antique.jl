# run `include("./developer/test.jl")`

dir = dirname(@__FILE__) * "/../"
cd(dir)
@show pwd()

using Pkg
Pkg.activate(dir)
Pkg.test()
