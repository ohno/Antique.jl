io = open("./result/InfinitePotentialWell3D.md", "w")


# <ψᵢ|ψⱼ> = ∫ψₙ*ψₙdx = δᵢⱼ


println(io, raw"""
#### Normalization & Orthogonality of $\psi_{n_x,n_y,n_z}(x,y,z)$

```math
\int_{0}^{L_x}\int_{0}^{L_y}\int_{0}^{L_z} \psi_{i_x,i_y,i_z}^\ast(x,y,z) \psi_{j_x,j_y,j_z}(x,y,z) ~\mathrm{d}x \mathrm{d}y\mathrm{d}z = \delta_{i_x,j_x}\delta_{i_y,j_y}\delta_{i_z,j_z}
```

```""")

@testset "IPW3D: <ψᵢ|ψⱼ> = ∫ψₙ*ψₙdx = δᵢⱼ" begin
  IPW3D = InfinitePotentialWell3D(L=[1.0,1.0,1.0], m=1.0, ℏ=1.0)
  @show IPW3D
  @show ψ(IPW3D,[0.5,0.5,0.5])
  for IPW3D in [
    InfinitePotentialWell3D(L=[1.0,1.0,1.0], m=1.0, ℏ=1.0)
    InfinitePotentialWell3D(L=[1.2,3.4,4.5], m=2.0, ℏ=3.0)
  ]
    println(io, "IPW3D = $IPW3D")
    println(io, "ix | iy | iz | jx | jy | jz |     analytical |      numerical ")
    println(io, "-- | -- | -- | -- | -- | -- | -------------- | -------------- ")
    for ix in 1:2
    for iy in 1:2
    for iz in 1:2
    for jx in 1:2
    for jy in 1:2
    for jz in 1:2
      Δr = 0.01
      numerical  = first(hcubature(r -> conj(ψ(IPW3D, r, n=[ix,iy,iz])) * ψ(IPW3D, r, n=[jx,jy,jz]), [Δr,Δr,Δr], IPW3D.L .- Δr, maxevals=1000))
      analytical = ((ix==jx && iy==jy && iz==jz) ? 1 : 0)
      acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-1) : isapprox(analytical, numerical, rtol=1e-1)
      @printf(io, "%2d | %2d | %2d | %2d | %2d | %2d | %14.9f | %14.9f %s\n", ix, iy, iz, jx, jy, jz, analytical, numerical, acceptance ? "✔" : "✗")
      @test acceptance
    end
    end
    end
    end
    end
    end
    println(io, "")
  end
end

println(io, """```\n""")


# <ψₙ|H|ψₙ>  = ∫ψₙ*Tψₙdx = Eₙ


println(io, raw"""
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

@testset "IPW3D: <ψₙ|H|ψₙ> = ∫ψₙ*Tψₙdx = Eₙ" begin
  IPW3D = InfinitePotentialWell3D(L=[1.0,1.0,1.0], m=1.0, ℏ=1.0)
  ∇²ψ(model,r;n=[1,1,1]) = sum(first(Zygote.diaghessian(x -> ψ(model,x,n=n), r)))
  ψTψ(model,r;n=[1,1,1]) = -model.ℏ^2/(2*model.m) * conj(ψ(model,r,n=n)) * ∇²ψ(model,r,n=n)
  @show IPW3D
  @show ψ(IPW3D,[0.5,0.5,0.5])
  @show ∇²ψ(IPW3D,[0.5,0.5,0.5])
  @show ψTψ(IPW3D,[0.5,0.5,0.5])
  for IPW3D in [
    InfinitePotentialWell3D(L=[1.0,1.0,1.0], m=1.0, ℏ=1.0)
    InfinitePotentialWell3D(L=[1.2,3.4,4.5], m=2.0, ℏ=3.0)
  ]
    println(io, "IPW3D = $IPW3D")
    println(io, " nx | ny | nz |     analytical |      numerical ")
    println(io, " -- | -- | -- | -------------- | -------------- ")
    for nx in [1,2]
    for ny in [1,2]
    for nz in [1,2]
      Δr = 0.01
      analytical = E(IPW3D, n=[nx,ny,nz])
      numerical  = first(hcubature(r -> ψTψ(IPW3D, r, n=[nx,ny,nz]), [Δr,Δr,Δr], IPW3D.L .- Δr, maxevals=1000))
      acceptance = iszero(analytical) ? isapprox(analytical, numerical, atol=1e-1) : isapprox(analytical, numerical, rtol=1e-1)
      @test acceptance
      @printf(io, " %2d | %2d | %2d | %14.9f | %14.9f %s\n", nx, ny, nz, numerical, analytical, acceptance ? "✔" :  "✗")
    end
    end
    end
    println(io, "")
  end
end

println(io, """```\n""")


close(io)