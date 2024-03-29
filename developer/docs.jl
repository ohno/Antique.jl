# run `include("./developer/docs.jl")`
# run(`julia --project=docs/ -e 'using Pkg; Pkg.activate("./"); cd("docs"); include("make.jl")'`)

dir = dirname(@__FILE__) * "/../"
@show dir
include(dir * "/docs/make.jl")
