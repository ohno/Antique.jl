module OldHydrogenAtom

    # Default
    Z = 1 # change here!
    ℏ = 1.0 # change here!
    Eₕ = 1.0 # change here!
    a₀ = 1.0 # change here!
    mₑ = 1.0 # change here!

    # Potential
    V(r; Z=Z, a₀=a₀) =  r<0 ? throw(DomainError(r, "The domain is 0≤r.")) : -Z/abs(r/a₀)

    # Energy
    E(; n=1, Z=Z, Eₕ=Eₕ) = -Z^2/(2*n^2) * Eₕ

    # Wave Function
    ψ(r, θ, φ; n=1, l=0, m=0, Z=Z, a₀=a₀) = R(r, n=n, l=l, Z=Z, a₀=a₀) * Y(θ, φ, l=l, m=m)

    # Radial Function
    function R(r; n=1, l=0, Z=Z, a₀=a₀)
        if !(0 ≤ r)
            throw(DomainError(r, "The domain is 0≤r."))
        end
        ρ = 2*Z*abs(r)/(n*a₀)
        N = -sqrt( factorial(n-l-1)/(2*n*factorial(n+l)) * (2*Z/(n*a₀))^3 )
        return N*ρ^l * exp(-ρ/2) * L(ρ, n=n+l, k=2*l+1)
    end

    # Associated Laguerre Polynomials
    L(x; n=0, k=0) = sum(m -> (-1)^(m+k) * factorial(n) // (factorial(m) * factorial(m+k) * factorial(n-m-k)) * x^m, 0:n-k)

    # Spherical Harmonics
    function Y(θ, φ; l=0, m=0)
        N = (-1)^((abs(m)+m)/2) * sqrt( (2*l+1)*factorial(l-Int(abs(m))) / (2*factorial(l+Int(abs(m)))) )
        return N * P(cos(θ),n=l,m=Int(abs(m))) * exp(im*m*φ) / sqrt(2*π)
    end

    # Associated Legendre Polynomials
    P(x; n=0, m=0) = (1//2)^n * (1-x^2)^(m//2) * sum(j -> (-1)^j * factorial(2*n-2*j) // (factorial(j) * factorial(n-j) * factorial(n-2*j-m)) * x^(n-2*j-m), 0:Int(floor((n-m)/2)))

end
