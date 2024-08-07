# Please run `include("./dev/logo.jl")` on RELP.

# https://github.com/JuliaLang/julia-logo-graphics
# https://products.aspose.app/imaging/conversion/svg-to-ico

dir = "$(@__DIR__)/../"
cd(dir)
@show pwd()
using Pkg
Pkg.activate(dir)
using Antique

function XY2path(X,Y)
    X = 252 / (maximum(X) - minimum(X)) * (X .- minimum(X)) .+ 50
    Y = 260 .- 445  .*  Y
    path = "M$(X[1]),$(Y[1]) "
    for i in 2:length(X)
        path *= "L$(X[i]),$(Y[i]) "
    end
    return path
end

HA = HydrogenAtom(Z=1, Eₕ=1.0, a₀=1.0, mₑ=1.0, ℏ=1.0)
paths = Dict()
X = 0.0:0.01:8
Y1 = [x^2 * Antique.R(HA,x,n=1,l=0)^2 for x in X]
Y2 = [x^2 * Antique.R(HA,x,n=2,l=0)^2 for x in X]
Y3 = [x^2 * Antique.R(HA,x,n=2,l=1)^2 for x in X]
Y4 = [x^2 * Antique.R(HA,x,n=3,l=0)^2 for x in X]
paths[1] = XY2path(X,Y1)
paths[2] = XY2path(X,Y2)
paths[3] = XY2path(X,Y3)
paths[4] = XY2path(X,Y4)

colors = ["#4063D8", "#CB3C33", "#389826", "#9558B2"]

as = 4
lw = 9

svg = """
<?xml version="1.0" encoding="UTF-8"?>
<svg
  version="1.1"
  xmlns="http://www.w3.org/2000/svg"
  xmlns:xlink="http://www.w3.org/1999/xlink"
  width="325pt"
  height="300pt"
  viewBox="0 0 325 300"
>

  <!-- Defs -->
  <defs>
    <marker id="arrow" markerWidth="$(as)" markerHeight="$(as)" refX="0" refY="$(as/2)" orient="auto">
      <polygon points="0,0 0,$(as) $(as),$(as/2)" fill="#24292E"/>
    </marker>
  </defs>

  <!-- Body --> 
  <path d="$(paths[4])" fill="none" stroke="$(colors[4])" stroke-width="$(lw)" stroke-linecap="round" />
  <path d="$(paths[3])" fill="none" stroke="$(colors[3])" stroke-width="$(lw)" stroke-linecap="round" />
  <path d="$(paths[2])" fill="none" stroke="$(colors[2])" stroke-width="$(lw)" stroke-linecap="round" />
  <path d="$(paths[1])" fill="none" stroke="$(colors[1])" stroke-width="$(lw)" stroke-linecap="round" />
  <path d="M10,260 L285,260," stroke="#24292E" stroke-width="$(lw)" stroke-linecap="round"  marker-end="url(#arrow)"/>
  <path d="M40,290 L40,40" stroke="#24292E" stroke-width="$(lw)" stroke-linecap="round"  marker-end="url(#arrow)"/>

</svg>
"""

HTML(svg) |> display

path = "./docs/src/assets/fig/logo.svg"
mkpath(dirname(path))
file = open(path, "w")
Base.write(file, svg)
close(file)
