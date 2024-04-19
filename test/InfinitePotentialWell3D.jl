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
  # for L in [0.1, 0.5, 1.0, 7.0]
  # for m in [0.1, 0.5, 1.0, 7.0]
  # for ℏ in [0.1, 0.5, 1.0, 7.0]
  for ix in 1:2
    for iy in 1:2
      for iz in 1:2
  for jx in 1:2
    for jy in 1:2
      for jz in 1:2
    analytical = ((ix==jx && iy==jy && iz==jz) ? 1 : 0)
    numerical  = quadgk(x -> 
                  quadgk(y ->
                    quadgk( z -> conj(ψ(IPW3D, x,y,z, nx=ix,ny=iy,nz=iz)) * ψ(IPW3D, x,y,z, nx=jx,ny=jy,nz=jz), 0.0, IPW3D.Lz, atol=10^-6)[1]
                    ,0.0, IPW3D.Ly, atol=10^-6)[1]
                    ,0.0, IPW3D.Lz, atol=10^-6)[1]
    acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-5) : isapprox(analytical, numerical, rtol=1e-5)
    @printf("%2d | %2d | %2d | %2d | %2d | %2d | %17.12f | %17.12f %s\n", ix, iy, iz, jx, jy, jz, analytical, numerical, acceptance ? "✔" : "✗")
    @test acceptance
      end
    end
  end
      end
    end
  end
end

println("""```
""")
