@doc raw"""
## Infinite Potential Well (Particle in a box)

### Schrödinger Equation
```math
  \psi(x) = E \psi(x)
```

### Hamiltonian
```math
   \frac{\hbar^2}{2m} \frac{\mathrm{d}^2}{\mathrm{d}x ^2} + V(x)
```

### Potential
`V(x; L=L)`
```math
  V(x) =
  \left\{
    \begin{array}{ll}
    \infty & x \lt 0, L \lt x \\
    0      & 0 \leq x \leq L
    \end{array}
  \right.
```

### Eigen Values
`E(; n=0, L=L, m=m, ℏ=ℏ)`
```math
  {\hbar^2 n^2 \pi^2}{2 m L^2}
```

### Eigen Functions
`ψ(x; n=0, L=L)`
```math
   \psi_n(x) = \sqrt{\frac{2}{L}} \sin \frac{n\pi x}{L}
```

### Proofs
- [Eigen Functions & Eigen Values](https://ja.wolframalpha.com/input?i2d=true&i=D%5B%5C%2840%29Sqrt%5BDivide%5B2%2Ca%5D%5Dsin%5C%2840%29Divide%5Bn%CF%80x%2Ca%5D%5C%2841%29%5C%2841%29%2C%7Bx%2C2%7D%5D)
- [Normalization](https://ja.wolframalpha.com/input?i=Integrate%5B%28%28Sqrt%5B2%2Fa%5Dsin%28%CF%80x%2Fa%29%29%29%5E2%2C+%7Bx%2C0%2Ca%7D%5D)
"""
module InfinitePotentialWell

  # Default
  L = 1.0
  m = 1.0 # 0.5
  ℏ = 1.0

  # Potential
  V(x; L=L) = 0<x<L ? 0 : Inf

  # Wave Function
  ψ(x; n=1, L=L) = sqrt(2/L) * sin(n*π*x/L)

  # Energy
  E(; n=1, L=L, m=m, ℏ=ℏ) = (ℏ^2*n^2*π^2) / (2*m*L^2)

end