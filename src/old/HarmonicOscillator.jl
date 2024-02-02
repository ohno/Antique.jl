module OldHarmonicOscillator

    # Default
    k = 1.0 # change here!
    m = 1.0 # change here!
    ℏ = 1.0 # change here!

    # Potential
    V(x; k=k, m=m, ω=sqrt(k/m)) = 1/2*k*x^2

    # Energy
    E(; n=0, k=k, m=m, ω=sqrt(k/m), ℏ=ℏ) = ℏ*ω*(n+1/2)

    # Wave Function
    function ψ(x; n=0, k=k, m=m, ω=sqrt(k/m), ℏ=ℏ)
        A = sqrt(1/(factorial(n)*2^n)*sqrt(m*ω/(π*ℏ)))
        ξ = sqrt(m*ω/ℏ) * x
        return A*H(ξ,n=n)*exp(-ξ^2/2)
    end

    # Hermite polynomials
    H(x; n=0) = factorial(n) * sum(i -> (-1)^i // (factorial(i)  * factorial(n-2*i)) * (2*x)^(n-2*i), 0:Int(floor(n/2)))

end
