

# <ψᵢ|ψⱼ> = ∫ψₙ*ψₙdx = δᵢⱼ


println(raw"""
#### Normalization & Orthogonality of $\psi_{n_x,n_y,n_z}(x,y,z)$

```math
\int_{0}^{L_x}\int_{0}^{L_y}\int_{0}^{L_z} \psi_{i_x,i_y,i_z}^\ast(x,y,z) \psi_{j_x,j_y,j_z}(x,y,z) ~\mathrm{d}x \mathrm{d}y\mathrm{d}z = \delta_{i_x,j_x}\delta_{i_y,j_y}\delta_{i_z,j_z}
```

```""")

@testset "<ψᵢ|ψⱼ> = ∫ψₙ*ψₙdx = δᵢⱼ" begin
  for IPW3D  in [
    InfinitePotentialWell3D(L=[1.0,1.0,1.0], m=1.0, ℏ=1.0)
    InfinitePotentialWell3D(L=[1.2,3.4,4.5], m=1.0, ℏ=1.0)
    InfinitePotentialWell3D(L=[1.2,3.4,4.5], m=2.0, ℏ=1.0)
    InfinitePotentialWell3D(L=[1.2,3.4,4.5], m=1.0, ℏ=2.0)
  ]
    @show IPW3D
    println("ix | iy | iz | jx | jy | jz |        analytical |         numerical ")
    println("-- | -- | -- | -- | -- | -- | ----------------- | ----------------- ")
    for ix in 1:2
    for iy in 1:2
    for iz in 1:2
    for jx in 1:2
    for jy in 1:2
    for jz in 1:2
      analytical = ((ix==jx && iy==jy && iz==jz) ? 1 : 0)
      numerical  = quadgk(x -> 
                   quadgk(y ->
                   quadgk(z ->
                     conj(ψ(IPW3D, x,y,z, n=[ix,iy,iz])) * ψ(IPW3D, x,y,z, n=[jx,jy,jz])
                   , 0.0, IPW3D.L[3], maxevals=10)[1]
                   , 0.0, IPW3D.L[2], maxevals=10)[1]
                   , 0.0, IPW3D.L[1], maxevals=10)[1]
      acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-5) : isapprox(analytical, numerical, rtol=1e-5)
      @printf("%2d | %2d | %2d | %2d | %2d | %2d | %17.12f | %17.12f %s\n", ix, iy, iz, jx, jy, jz, analytical, numerical, acceptance ? "✔" : "✗")
      @test acceptance
    end
    end
    end
    end
    end
    end
    println()
  end
end

println("""```
""")



# <ψₙ|H|ψₙ>  = ∫ψₙ*Tψₙdx = Eₙ


println(raw"""
#### Eigenvalues

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

ψTψ(IPW3D, x, y, z; n=[1,1,1], Δx=0.01, Δy=0.01, Δz=0.01) = -IPW3D.ℏ^2/(2*IPW3D.m) * conj(ψ(IPW3D,x,y,z,n=n)) * (
  ( ψ(IPW3D,x+Δx,y,z,n=n) -2*ψ(IPW3D,x,y,z,n=n) + ψ(IPW3D,x-Δx,y,z,n=n) ) / Δx^2 +
  ( ψ(IPW3D,x,y+Δy,z,n=n) -2*ψ(IPW3D,x,y,z,n=n) + ψ(IPW3D,x,y-Δy,z,n=n) ) / Δy^2 +
  ( ψ(IPW3D,x,y,z+Δz,n=n) -2*ψ(IPW3D,x,y,z,n=n) + ψ(IPW3D,x,y,z-Δz,n=n) ) / Δz^2
)

@testset "<ψₙ|H|ψₙ>  = ∫ψₙ*Tψₙdx = Eₙ" begin
  for IPW3D  in [
    InfinitePotentialWell3D(L=[1.0,1.0,1.0], m=1.0, ℏ=1.0)
    InfinitePotentialWell3D(L=[1.2,3.4,4.5], m=1.0, ℏ=1.0)
    InfinitePotentialWell3D(L=[1.2,3.4,4.5], m=2.0, ℏ=1.0)
    InfinitePotentialWell3D(L=[1.2,3.4,4.5], m=1.0, ℏ=2.0)
  ]
    @show IPW3D
    println(" nx | ny | nz |        analytical |         numerical ")
    println(" -- | -- | -- | ----------------- | ----------------- ")
    for nx in [1,2]
    for ny in [1,2]
    for nz in [1,2]
      analytical = E(IPW3D,n=[nx,ny,nz])
      numerical  = quadgk(x ->
                   quadgk(y ->
                   quadgk(z ->
                     ψTψ(IPW3D, x, y, z, n=[nx,ny,nz], Δx=IPW3D.L[1]*0.0001, Δy=IPW3D.L[2]*0.0001, Δz=IPW3D.L[3]*0.0001)
                   , 0, IPW3D.L[3], maxevals=5)[1]
                   , 0, IPW3D.L[2], maxevals=5)[1]
                   , 0, IPW3D.L[1], maxevals=5)[1]
      acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-5) : isapprox(analytical, numerical, rtol=1e-5)
      @test acceptance
      @printf(" %2d | %2d | %2d | %17.12f | %17.12f %s\n", nx, ny, nz, numerical, analytical, acceptance ? "✔" :  "✗")
    end
    end
    end
    println()
  end
end

println("""```""")