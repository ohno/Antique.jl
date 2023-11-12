@doc raw"""
## Harmonic Oscillator

### Schrödinger Equation
```math
    \hat{H}\psi(x) = E \psi(x)
```

### Hamiltonian
```math
    \begin{aligned}
        \hat{H}
        &= - \frac{\hbar^2}{2m} \frac{\mathrm{d}^2}{\mathrm{d}x ^2} + V(x) \\
        &= - \frac{1}{2} \hbar\omega \frac{\mathrm{d}^2}{\mathrm{d}\xi^2} + V(x)
    \end{aligned}
```

### Potential
`V(x; k=k, m=m, ω=sqrt(k/m))`
```math
    V(x)
    = \frac{1}{2} k x^2
    = \frac{1}{2} m \omega^2 x^2
    = \frac{1}{2} \hbar \omega \xi^2
```

### Eigen Values
`E(n; k, m, ω=sqrt(k/m), ℏ=ℏ)`
```math
    E_n = \hbar \omega \left( n + \frac{1}{2} \right)
```

where
```math
    \omega = \sqrt{k/m}.
```

### Eigen Functions
`ψ(n, x; k=k, m=m, ω=sqrt(k/m), ℏ=ℏ)`
```math
    \psi_n(x) = A_n H_n(\xi) \exp{\left( -\frac{\xi^2}{2} \right)}
```
where
```math
    \omega = \sqrt{k/m},
    \xi = \sqrt{\frac{m\omega}{\hbar}}x,
    A_n = \sqrt{\frac{1}{n! 2^n} \sqrt{\frac{m\omega}{\pi\hbar}}},
    H_n(x) = (-1)^n \mathrm{e}^{x^2} \frac{\mathrm{d}^n}{\mathrm{d}x^n} \mathrm{e}^{-x^2}.
```

### Hermite Polynomials
`H(x; n=0)`

Rodrigues' formula & closed-form:
```math
    \begin{aligned}
        H_{n}(x)
        &:= (-1)^n \mathrm{e}^{x^2} \frac{\mathrm{d}^n}{\mathrm{d}x^n} \mathrm{e}^{-x^2} \\
        &= n! \sum_{m=0}^{\lfloor n/2 \rfloor} \frac{(-1)^m}{m! (n-2m)!}(2 x)^{n-2m}.
    \end{aligned}
```
Examples:
```math
    \begin{aligned}
        H_{0}(x)  &= 1, \\
        H_{1}(x)  &= 2 x, \\
        H_{2}(x)  &= -2 + 4 x^{2}, \\
        H_{3}(x)  &= -12 x + 8 x^{3}, \\
        H_{4}(x)  &= 12 - 48 x^{2} + 16 x^{4}, \\
        H_{5}(x)  &= 120 x - 160 x^{3} + 32 x^{5}, \\
        H_{6}(x)  &= -120 + 720 x^{2} - 480 x^{4} + 64 x^{6}, \\
        H_{7}(x)  &= -1680 x + 3360 x^{3} - 1344 x^{5} + 128 x^{7}, \\
        H_{8}(x)  &= 1680 - 13440 x^{2} + 13440 x^{4} - 3584 x^{6} + 256 x^{8}, \\
        H_{9}(x)  &= 30240 x - 80640 x^{3} + 48384 x^{5} - 9216 x^{7} + 512 x^{9}, \\
        &\vdots
    \end{aligned}
```

### Reference
- [DLMF 18.5.18](https://dlmf.nist.gov/18.5#E18)
- [cpprefjp](https://cpprefjp.github.io/reference/cmath/hermite.html)
- The Digital Library of Mathematical Functions (DLMF)                                                    [18.3 Table1](https://dlmf.nist.gov/18.3#T1), [18.5 Table1](https://dlmf.nist.gov/18.5#T1), [18.5.13](https://dlmf.nist.gov/18.5#E13), [18.5.18](https://dlmf.nist.gov/18.5#E18)
- L. D. Landau, E. M. Lifshitz, Quantum Mechanics (Pergamon Press, 1965)                                  [p.595 (a.4), (a.6)](https://archive.org/details/ost-physics-landaulifshitz-quantummechanics/page/n607/mode/2up)
- L. I. Schiff, Quantum Mechanics (McGraw-Hill Book Company, 1968)                                        [p.71 (13.12)](https://archive.org/details/ost-physics-schiff-quantummechanics/page/n87/mode/2up)
- A. Messiah, Quanfum Mechanics (Dover Publications, 1999)                                                [p.491 (B.59)](https://archive.org/details/quantummechanics0000mess/page/491/mode/1up)
- W. Greiner, Quantum Mechanics: An Introduction Third Edition (Springer, 1994)                           [p.152 (7.22)](https://archive.org/details/quantummechanics0001grei_u4x0/page/152/mode/1up)
- D. J. Griffiths, Introduction to Quantum Mechanics (Prentice Hall, 1995)                                [p.41 Table 2.1](https://archive.org/details/griffiths-introduction-to-quantum-mechanics/page/41/mode/1up), [p.43 (2.70)](https://archive.org/details/griffiths-introduction-to-quantum-mechanics/page/43/mode/1up)
- D. A. McQuarrie, J. D. Simon, Physical Chemistry: A Molecular Approach (University Science Books, 1997) [p.170 Table 5.2](https://archive.org/details/McQuarrieSimonPhysicalChemistrySolutions/McQuarrie_Simon_Physical_Chemistry1997/page/n193/mode/1up)
- P. W. Atkins, J. De Paula, Atkins' Physical Chemistry, 8th edition (W. H. Freeman, 2008)                [p.293 Table 9.1](https://archive.org/details/atkinsphysicalch00pwat/page/292/mode/2up)
- J. J. Sakurai, J. Napolitano, Modern Quantum Mechanics Third Edition (Cambridge University Press, 2021) [p.524 (B.29)](https://doi.org/10.1017/9781108587280)
"""
module HarmonicOscillator

    # Default
    k = 1.0
    m = 1.0 # 3.0
    ℏ = 1.0

    # Packages
    # using SpecialPolynomials

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
    H(x; n=0) = factorial(n) * sum(m -> (-1)^m // (factorial(m)  * factorial(n-2*m)) * (2*x)^(n-2*m), 0:Int(floor(n/2)))

end
