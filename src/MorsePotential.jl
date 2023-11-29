module MorsePotential

    # Packages
    using SpecialFunctions

    # Default
    # F. M. Fernández, J. Garcia, ChemistrySelect, 6, 9527−9534(2021)
    # https://doi.org/10.1002/slct.202102509
    # https://physics.nist.gov/cgi-bin/cuu/Value?mpsme
    rₑ =  1.997193319969992120068298141276 # change here!
    Vₑ = -0.602634619106539878727562156289 # change here!
    Dₑ = - 0.5 - Vₑ # change here!
    k = 2*((-1.1026342144949464615+1/2.00) - Vₑ) / (2.00 - rₑ)^2 # change here!
    µ = 1/(1/1836.15267343 + 1/1836.15267343) # change here!
    ℏ = 1.0 # change here!

    # Potential
    V(r; rₑ=rₑ, Dₑ=Dₑ, k=k, a=sqrt(k/(2*Dₑ))) = r<0 ? throw(DomainError(r, "r=$r is out of the domain (0≦r)")) : Dₑ*( exp(-2*a*(r-rₑ)) -2*exp(-a*(r-rₑ)) )

    # Energy
    E(; n=0, rₑ=rₑ, Dₑ=Dₑ, k=k, a=sqrt(k/(2*Dₑ)), µ=µ, ω=sqrt(k/µ), χ=ℏ*ω/(4*Dₑ), ℏ=ℏ) = - Dₑ + ℏ*ω*(n+1/2) - χ*ℏ*ω*(n+1/2)^2

    # Maximum of n
    nₘₐₓ(; Dₑ=Dₑ, k=k, µ=µ, ω=sqrt(k/µ)) = Int(floor((2*Dₑ - ω)/ω))

    # Wave Function
    function ψ(r; n=0, rₑ=rₑ, Dₑ=Dₑ, k=k, a=sqrt(k/(2*Dₑ)), µ=µ, ω=sqrt(k/µ), χ=ℏ*ω/(4*Dₑ), ℏ=ℏ)
        if r<0
            throw(DomainError(r, "r=$r is out of the domain (0≦r)"))
        end
        λ = sqrt(2*µ*Dₑ) / (a*ℏ)
        ξ = 2*λ*exp(-a*(r-rₑ))
        s  = 2*λ - 2*n - 1
        N  = sqrt(factorial(n) * s * a / Γ(s+n+1))
        return N * ξ^(s/2) * exp(-ξ/2) * L(ξ,n=n,α=s)
    end

    # Gamma Function
    Γ(z) = try; gamma(z); catch; Inf; end

    # Generalized Laguerre Polynomials
    L(x; n=0, α=0) = sum(k -> (-1)^(k) * (gamma(α+n+1) / (gamma(α+1+k)*gamma(n-k+1))) * x^k / factorial(k), 0:n)
    Lαint(x; n=0, α::Int=0) = sum(k -> (-1)^(k) * (Int(gamma(α+n+1)) // Int((gamma(α+1+k)*gamma(n-k+1)))) * x^k // factorial(k), 0:n) # α::Int

end