# Please run `include("./dev/revice.jl")` on RELP.

using Pkg
# Pkg.instantiate()

try
  using Revise
catch
  Pkg.add("Revise")
  using Revise
end

dir = dirname(@__FILE__) * "/../"
cd(dir)
@show pwd()
Pkg.activate(dir)
using Antique
