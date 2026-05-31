#### Hermite Polynomials $H_n(x)$

```math
  \begin{aligned}
    H_{n}(x)
    &:= (-1)^n \mathrm{e}^{x^2} \frac{\mathrm{d}^n}{\mathrm{d}x^n} \mathrm{e}^{-x^2} \\
    &= n! \sum_{m=0}^{\lfloor n/2 \rfloor} \frac{(-1)^m}{m! (n-2m)!}(2 x)^{n-2m}.
  \end{aligned}
```

#### Normalization & Orthogonality of $H_n(x)$

```math
\int_{-\infty}^\infty H_j(x) H_i(x) \mathrm{e}^{-x^2} \mathrm{d}x = \sqrt{\pi} 2^j j! \delta_{ij}
```

```
 i |  j |     analytical |      numerical 
-- | -- | -------------- | -------------- 
```

#### Normalization & Orthogonality of $\psi_n(x)$

```math
\int \psi_i^\ast(x) \psi_j(x) \mathrm{d}x = \delta_{ij}
```

```
 i |  j |     analytical |      numerical 
-- | -- | -------------- | -------------- 
```

#### Virial Theorem

The virial theorem $\langle T \rangle = \langle V \rangle$ and the definition of Hamiltonian $\langle H \rangle = \langle T \rangle + \langle V \rangle$ derive $\langle H \rangle = 2 \langle V \rangle = 2 \langle T \rangle$.

```math
2 \int \psi_n^\ast(x) V(x) \psi_n(x) \mathrm{d}x = E_n
```

```
  k |  n |     analytical |      numerical 
--- | -- | -------------- | -------------- 
```

#### Eigenvalues

```math
  \begin{aligned}
    E_n
    &=      \int \psi^\ast_n(x) \hat{H} \psi_n(x) \mathrm{d}x \\
    &=      \int \psi^\ast_n(x) \left[ \hat{V} + \hat{T} \right] \psi(x) \mathrm{d}x \\
    &=      \int \psi^\ast_n(x) \left[ V(x) - \frac{\hbar^2}{2m} \frac{\mathrm{d}^{2}}{\mathrm{d} x^{2}} \right] \psi(x) \mathrm{d}x \\
    &\simeq \int \psi^\ast_n(x) \left[ V(x)\psi(x) -\frac{\hbar^2}{2m} \frac{\psi(x+\Delta x) - 2\psi(x) + \psi(x-\Delta x)}{\Delta x^{2}} \right] \mathrm{d}x.
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

```
  k |  n |     analytical |      numerical 
--- | -- | -------------- | -------------- 
```

