@doc raw"""
## Morse Potential

### Schrodinger Equation
```math
    \hat{H}\psi(r) = E \psi(r)
```

### Hamiltonian
```math
    \hat{H} = - \frac{\hbar^2}{2\mu} \frac{\mathrm{d}^2}{\mathrm{d}r ^2} + V(r)
```

### Potential
`V(r; rₑ=rₑ, Dₑ=Dₑ, k=k, a=sqrt(k/(2*Dₑ)))`
```math
    V(r) = D_\mathrm{e} \left( \mathrm{e}^{-2a(r-r_e)} - 2\mathrm{e}^{-a(r-r_e)} \right)
```

### Eigen Values
`E(n; rₑ=rₑ, Dₑ=Dₑ, k=k, a=sqrt(k/(2*Dₑ)), µ=µ, ω=sqrt(k/µ), χ=ℏ*ω/(4*Dₑ), ℏ=ℏ)`
```math
    E_n = - D_\mathrm{e} + \hbar \omega \left( n + \frac{1}{2} \right) - \chi \hbar \omega \left( n + \frac{1}{2} \right)^2
```

### Eigen Functions
`ψ(n, r; rₑ=rₑ, Dₑ=Dₑ, k=k, a=sqrt(k/(2*Dₑ)), µ=µ, ω=sqrt(k/µ), χ=ℏ*ω/(4*Dₑ), ℏ=ℏ)`
```math
    \psi_n(r) = N_n z^{\lambda-n-1/2} \mathrm{e}^{-z/2} L_n^{(2\lambda-2n-1)}(\xi)
```

where, 
```math
    \xi := 2\lambda\mathrm{e}^{-a(r-r_e)},
    \omega := \sqrt{k/µ},
    k := 2D_\mathrm{e}a^2,
    \lambda := \frac{\sqrt{2mD_\mathrm{e}}}{a\hbar},
    \chi := \frac{\hbar\omega}{4D_\mathrm{e}},
    N_n := \sqrt{\frac{n!(2\lambda-2n-1)a}{\Gamma(2\lambda-n)}},
    L_n^{(\alpha)}(x) := \frac{x^{-\alpha} \mathrm{e}^x}{n !} \frac{\mathrm{d}^n}{\mathrm{d} x^n}\left(\mathrm{e}^{-x} x^{n+\alpha}\right).
```

### Generalized Laguerre Polynomials
`L(x; n=0, α=0)`

Rodrigues' formula & closed-form:
```math
    \begin{aligned}
        L_n^{(\alpha)}(x)
        &= \frac{x^{-\alpha}e^x}{n!} \frac{d^n}{dx^n}\left(x^{n+\alpha}e^{-x}\right) \\
        &= \sum_{k=0}^n(-1)^k \left(\begin{array}{l} n+\alpha \\ n-k \end{array}\right) \frac{x^k}{k !} \\
        &= \sum_{k=0}^n(-1)^k \frac{\Gamma(\alpha+n+1)}{\Gamma(\alpha+k+1)\Gamma(n-k+1)} \frac{x^k}{k !}.
    \end{aligned}
```
Examples:
```math
    \begin{aligned}
        L_0^{(0)}(x) &= 1, \\
        L_1^{(0)}(x) &= 1 - x, \\
        L_1^{(1)}(x) &= 2 - x, \\
        L_2^{(0)}(x) &= 1 - 2 x + 1/2 x^{2}, \\
        L_2^{(1)}(x) &= 3 - 3 x + 1/2 x^{2}, \\
        L_2^{(2)}(x) &= 6 - 4 x + 1/2 x^{2}, \\
        L_3^{(0)}(x) &= 1 - 3 x + 3/2 x^{2} - 1/6 x^{3}, \\
        L_3^{(1)}(x) &= 4 - 6 x + 2 x^{2} - 1/6 x^{3}, \\
        L_3^{(2)}(x) &= 10 - 10 x + 5/2 x^{2} - 1/6 x^{3}, \\
        L_3^{(3)}(x) &= 20 - 15 x + 3 x^{2} - 1/6 x^{3}, \\
        L_4^{(0)}(x) &= 1 - 4 x + 3 x^{2} - 2/3 x^{3} + 1/24 x^{4}, \\
        L_4^{(1)}(x) &= 5 - 10 x + 5 x^{2} - 5/6 x^{3} + 1/24 x^{4}, \\
        L_4^{(2)}(x) &= 15 - 20 x + 15/2 x^{2} - 1 x^{3} + 1/24 x^{4}, \\
        L_4^{(3)}(x) &= 35 - 35 x + 21/2 x^{2} - 7/6 x^{3} + 1/24 x^{4}, \\
        L_4^{(4)}(x) &= 70 - 56 x + 14 x^{2} - 4/3 x^{3} + 1/24 x^{4}, \\
        \vdots
    \end{aligned}
```

### References
- [P. M. Morse, Phys. Rev. 34, 57 (1929)](https://doi.org/10.1103/PhysRev.34.57)
- [J. P. Dahl, M. Springborg, J. Chem. Phys. 88, 4535 (1988). (62), (63)](https://doi.org/10.1063/1.453761)
- [W. K. Shao, Y. He, J. Pan, J. Nonlinear Sci. Appl., 9, 5, 3388 (2016). (1.6)](http://dx.doi.org/10.22436/jnsa.009.05.124) 
- The Digital Library of Mathematical Functions (DLMF) [18.3 Table1](https://dlmf.nist.gov/18.3#T1), [18.5 Table1](https://dlmf.nist.gov/18.5#T1), [18.5.12](https://dlmf.nist.gov/18.5#E12)
"""
module MorsePotential

    # Default
    # F. M. Fernández, J. Garcia, ChemistrySelect, 6, 9527−9534(2021)
    # https://doi.org/10.1002/slct.202102509
    # https://physics.nist.gov/cgi-bin/cuu/Value?mpsme
    rₑ =  1.997193319969992120068298141276
    Vₑ = -0.602634619106539878727562156289
    Dₑ = - 0.5 - Vₑ    
    k = 2*((-1.1026342144949464615+1/2.00) - Vₑ) / (2.00 - rₑ)^2
    µ = 1/(1/1836.15267343 + 1/1836.15267343)
    ℏ = 1.0

    # Packages
    using SpecialFunctions

    # Potential
    V(r; rₑ=rₑ, Dₑ=Dₑ, k=k, a=sqrt(k/(2*Dₑ))) = Dₑ*( exp(-2*a*(r-rₑ)) -2*exp(-a*(r-rₑ)) )

    # Energy
    E(; n=0, rₑ=rₑ, Dₑ=Dₑ, k=k, a=sqrt(k/(2*Dₑ)), µ=µ, ω=sqrt(k/µ), χ=ℏ*ω/(4*Dₑ), ℏ=ℏ) = - Dₑ + ℏ*ω*(n+1/2) - χ*ℏ*ω*(n+1/2)^2

    # Maximum of n
    nₘₐₓ(; Dₑ=Dₑ, k=k, µ=µ, ω=sqrt(k/µ)) = Int(floor((2*Dₑ - ω)/ω))

    # Wave Function
    function ψ(r; n=0, rₑ=rₑ, Dₑ=Dₑ, k=k, a=sqrt(k/(2*Dₑ)), µ=µ, ω=sqrt(k/µ), χ=ℏ*ω/(4*Dₑ), ℏ=ℏ)
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