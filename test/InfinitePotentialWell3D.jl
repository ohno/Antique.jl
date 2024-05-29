IPW3D = InfinitePotentialWell3D(Lx=1.0, Ly=1.1, Lz=1.2, m=1.0, ℏ=1.0)


# <ψᵢ|ψⱼ> = ∫ψₙ*ψₙdx = δᵢⱼ


println(raw"""
#### Normalization & Orthogonality of $\psi_{n_x,n_y,n_z}(x,y,z)$

```math
\int_{0}^{L_x}\int_{0}^{L_y}\int_{0}^{L_z} \psi_{i_x,i_y,i_z}^\ast(x,y,z) \psi_{j_x,j_y,j_z}(x,y,z) ~\mathrm{d}x \mathrm{d}y\mathrm{d}z = \delta_{i_x,j_x}\delta_{i_y,j_y}\delta_{i_z,j_z}
```

```""")

@testset "<ψᵢ|ψⱼ> = ∫ψₙ*ψₙdx = δᵢⱼ" begin
  println("ix | iy | iz | jx | jy | jz |        analytical |         numerical ")
  println("-- | -- | -- | -- | -- | -- | ----------------- | ----------------- ")
  # for ix in 1:2
  # for iy in 1:2
  # for iz in 1:2
  # for jx in 1:2
  # for jy in 1:2
  # for jz in 1:2
  #   analytical = ((ix==jx && iy==jy && iz==jz) ? 1 : 0)
  #   numerical  = quadgk(x -> 
  #                quadgk(y ->
  #                quadgk(z ->
  #                  conj(ψ(IPW3D, x,y,z, nx=ix, ny=iy, nz=iz)) * ψ(IPW3D, x,y,z, nx=jx, ny=jy, nz=jz)
  #                , 0.0, IPW3D.Lz, maxevals=10)[1]
  #                , 0.0, IPW3D.Ly, maxevals=10)[1]
  #                , 0.0, IPW3D.Lz, maxevals=10)[1]
  #   acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-2) : isapprox(analytical, numerical, rtol=1e-2)
  #   @printf("%2d | %2d | %2d | %2d | %2d | %2d | %17.12f | %17.12f %s\n", ix, iy, iz, jx, jy, jz, analytical, numerical, acceptance ? "✔" : "✗")
  #   @test acceptance
  # end
  # end
  # end
  # end
  # end
  # end
end

println("""```
""")



# <ψₙ|H|ψₙ>  = ∫ψₙ*Tψₙdx = Eₙ


println(raw"""
#### Eigen Values

```math
  \begin{aligned}
    E_n
    &=      \int_0^L \psi^\ast_n(x) \hat{H} \psi_n(x) ~\mathrm{d}x \\
    &=      \int_0^L \psi^\ast_n(x) \left[ \hat{V} + \hat{T} \right] \psi(x) ~\mathrm{d}x \\
    &=      \int_0^L \psi^\ast_n(x) \left[ 0 - \frac{\hbar^2}{2m} \frac{\mathrm{d}^{2}}{\mathrm{d} x^{2}} \right] \psi(x) ~\mathrm{d}x \\
    &\simeq \int_0^L \psi^\ast_n(x) \left[ -\frac{\hbar^2}{2m} \frac{\psi(x+\Delta x) - 2\psi(x) + \psi(x-\Delta x)}{\Delta x^{2}} \right] ~\mathrm{d}x.
  \end{aligned}
```

Where, the difference formula for the 2nd-order derivative:

```math
\begin{aligned}
  % 2\psi(x)
  % + \frac{\mathrm{d}^{2} \psi(x)}{\mathrm{d} x^{2}} \Delta x^{2}
  % + O\left(\Delta x^{4}\right)
  % &=
  % \psi(x+\Delta x)
  % + \psi(x-\Delta x)
  % \\
  % \frac{\mathrm{d}^{2} \psi(x)}{\mathrm{d} x^{2}} \Delta x^{2}
  % &=
  % \psi(x+\Delta x)
  % - 2\psi(x)
  % + \psi(x-\Delta x)
  % - O\left(\Delta x^{4}\right)
  % \\
  % \frac{\mathrm{d}^{2} \psi(x)}{\mathrm{d} x^{2}}
  % &=
  % \frac{\psi(x+\Delta x) - 2\psi(x) + \psi(x-\Delta x)}{\Delta x^{2}}
  % - \frac{O\left(\Delta x^{4}\right)}{\Delta x^{2}}
  % \\
  \frac{\mathrm{d}^{2} \psi(x)}{\mathrm{d} x^{2}}
  &=
  \frac{\psi(x+\Delta x) - 2\psi(x) + \psi(x-\Delta x)}{\Delta x^{2}}
  + O\left(\Delta x^{2}\right)
\end{aligned}
```

are given by the sum of 2 Taylor series:

```math
\begin{aligned}
\psi(x+\Delta x)
&= \psi(x)
+ \frac{\mathrm{d} \psi(x)}{\mathrm{d} x} \Delta x
+ \frac{1}{2!} \frac{\mathrm{d}^{2} \psi(x)}{\mathrm{d} x^{2}} \Delta x^{2}
+ \frac{1}{3!} \frac{\mathrm{d}^{3} \psi(x)}{\mathrm{d} x^{3}} \Delta x^{3}
+ O\left(\Delta x^{4}\right),
\\
\psi(x-\Delta x)
&= \psi(x)
- \frac{\mathrm{d} \psi(x)}{\mathrm{d} x} \Delta x
+ \frac{1}{2!} \frac{\mathrm{d}^{2} \psi(x)}{\mathrm{d} x^{2}} \Delta x^{2}
- \frac{1}{3!} \frac{\mathrm{d}^{3} \psi(x)}{\mathrm{d} x^{3}} \Delta x^{3}
+ O\left(\Delta x^{4}\right).
\end{aligned}
```
```""")

ψTψ(IPW3D, x,y,z; nx=0,ny=0,nz=0, Δx=0.01,Δy=0.01,Δz=0.01) = -IPW3D.ℏ^2/(2*IPW3D.m) * conj(ψ(IPW3D,x,y,z,nx=nx,ny=ny,nz=nz)) * (
  ( ψ(IPW3D,x+Δx,y,z,nx=nx,ny=ny,nz=nz) -2*ψ(IPW3D,x,y,z,nx=nx,ny=ny,nz=nz) + ψ(IPW3D,x-Δx,y,z,nx=nx,ny=ny,nz=nz) ) / Δx^2 +
  ( ψ(IPW3D,x,y+Δy,z,nx=nx,ny=ny,nz=nz) -2*ψ(IPW3D,x,y,z,nx=nx,ny=ny,nz=nz) + ψ(IPW3D,x,y-Δy,z,nx=nx,ny=ny,nz=nz) ) / Δy^2 +
  ( ψ(IPW3D,x,y,z+Δz,nx=nx,ny=ny,nz=nz) -2*ψ(IPW3D,x,y,z,nx=nx,ny=ny,nz=nz) + ψ(IPW3D,x,y,z-Δz,nx=nx,ny=ny,nz=nz) ) / Δz^2
)

@testset "<ψₙ|H|ψₙ>  = ∫ψₙ*Tψₙdx = Eₙ" begin
  println(" nx |  ny |  nz |        analytical |         numerical ")
  println(" -- | --- | --- | ----------------- | ----------------- ")
  # for nx in [1,2]
  # for ny in [1,2]
  # for nz in [1,2]
  #   IPW3D = InfinitePotentialWell3D(Lx=1.0,Ly=2.0,Lz=3.0)
  #   analytical = E(IPW3D,nx=nx,ny=ny,nz=nz)
  #   numerical  = quadgk(x ->
  #                quadgk(y ->
  #                quadgk(z ->
  #                  ψTψ(IPW3D, x, y, z, nx=nx, ny=ny, nz=nz, Δx=IPW3D.Lx*0.0001, Δy=IPW3D.Ly*0.0001, Δz=IPW3D.Lz*0.0001)
  #                , 0, IPW3D.Lz, maxevals=5)[1]
  #                , 0, IPW3D.Ly, maxevals=5)[1]
  #                , 0, IPW3D.Lz, maxevals=5)[1]
  #   acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-1) : isapprox(analytical, numerical, rtol=1e-1)
  #   @test acceptance
  #   @printf(" %2d | %3d | %3d | %17.12f | %17.12f %s\n", nx, ny, nz, numerical, analytical, acceptance ? "✔" :  "✗")
  # end
  # end
  # end
end

println("""```""")