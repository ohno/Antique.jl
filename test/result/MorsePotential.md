#### Generalized Laguerre Polynomials $L_n^{(\alpha)}(x)$

```math
  \begin{aligned}
    L_n^{(\alpha)}(x)
    &= \frac{x^{-\alpha}e^x}{n!} \frac{d^n}{dx^n}\left(x^{n+\alpha}e^{-x}\right) \\
    &= \sum_{k=0}^n(-1)^k \frac{\Gamma(\alpha+n+1)}{\Gamma(\alpha+k+1)\Gamma(n-k+1)} \frac{x^k}{k !}.
  \end{aligned}
```

#### Normalization & Orthogonality of $L_n^{(\alpha)}(x)$

```math
\int_0^\infty L_i^{(\alpha)}(x) L_j^{(\alpha)}(x) x^\alpha \mathrm{e}^{-x} \mathrm{d}x = \frac{\Gamma(n+\alpha+1)}{n!} \delta_{ij}
```

```
   α |  i |  j |     analytical |      numerical 
---- | -- | -- | -------------- | -------------- 
```

#### Normalization & Orthogonality of $\psi_n(r)$

```math
\int_0^\infty \psi_i^\ast(r) \psi_j(r) \mathrm{d}r = \delta_{ij}
```

```
 i |  j |     analytical |      numerical 
-- | -- | -------------- | -------------- 
```

#### Eigenvalues

```math
  \begin{aligned}
    E_n
    &=      \int \psi^\ast_n(r) \hat{H} \psi_n(r) \mathrm{d}r \\
    &=      \int \psi^\ast_n(r) \left[ \hat{V} + \hat{T} \right] \psi(r) \mathrm{d}r \\
    &=      \int \psi^\ast_n(r) \left[ V(r) - \frac{\hbar^2}{2m} \frac{\mathrm{d}^{2}}{\mathrm{d} r^{2}} \right] \psi(r) \mathrm{d}r \\
    &\simeq \int \psi^\ast_n(r) \left[ V(r)\psi(r) -\frac{\hbar^2}{2m} \frac{\psi(r+\Delta r) - 2\psi(r) + \psi(r-\Delta r)}{\Delta r^{2}} \right] \mathrm{d}r.
  \end{aligned}
```

Where, the difference formula for the 2nd-order derivative:

```math
\begin{aligned}
  % 2\psi(r)
  % + \frac{\mathrm{d}^{2} \psi(r)}{\mathrm{d} r^{2}} \Delta r^{2}
  % + O\left(\Delta r^{4}\right)
  % &=
  % \psi(r+\Delta r)
  % + \psi(r-\Delta r)
  % \\
  % \frac{\mathrm{d}^{2} \psi(r)}{\mathrm{d} r^{2}} \Delta r^{2}
  % &=
  % \psi(r+\Delta r)
  % - 2\psi(r)
  % + \psi(r-\Delta r)
  % - O\left(\Delta r^{4}\right)
  % \\
  % \frac{\mathrm{d}^{2} \psi(r)}{\mathrm{d} r^{2}}
  % &=
  % \frac{\psi(r+\Delta r) - 2\psi(r) + \psi(r-\Delta r)}{\Delta r^{2}}
  % - \frac{O\left(\Delta r^{4}\right)}{\Delta r^{2}}
  % \\
  \frac{\mathrm{d}^{2} \psi(r)}{\mathrm{d} r^{2}}
  &=
  \frac{\psi(r+\Delta r) - 2\psi(r) + \psi(r-\Delta r)}{\Delta r^{2}}
  + O\left(\Delta r^{2}\right)
\end{aligned}
```

are given by the sum of 2 Taylor series:

```math
\begin{aligned}
\psi(r+\Delta r)
&= \psi(r)
+ \frac{\mathrm{d} \psi(r)}{\mathrm{d} r} \Delta r
+ \frac{1}{2!} \frac{\mathrm{d}^{2} \psi(r)}{\mathrm{d} r^{2}} \Delta r^{2}
+ \frac{1}{3!} \frac{\mathrm{d}^{3} \psi(r)}{\mathrm{d} r^{3}} \Delta r^{3}
+ O\left(\Delta r^{4}\right),
\\
\psi(r-\Delta r)
&= \psi(r)
- \frac{\mathrm{d} \psi(r)}{\mathrm{d} r} \Delta r
+ \frac{1}{2!} \frac{\mathrm{d}^{2} \psi(r)}{\mathrm{d} r^{2}} \Delta r^{2}
- \frac{1}{3!} \frac{\mathrm{d}^{3} \psi(r)}{\mathrm{d} r^{3}} \Delta r^{3}
+ O\left(\Delta r^{4}\right).
\end{aligned}
```

```
  k |  n |     analytical |      numerical 
--- | -- | -------------- | -------------- 
```

#### Recurrence Relation between $E_{n+1}$ and $E_n$

```math
\begin{equation}
\left\{ \,
  \begin{aligned}
    0 < \Delta E && 0 \leq n \leq n_\mathrm{max} \\
    \Delta E < 0 && \mathrm{otherwise}
  \end{aligned}
\right.
\end{equation}
```

```math
\Delta E =  E_{n+1} - E_n
```

```math
n_\mathrm{max} = \left\lfloor\frac{2 D_{\mathrm{e}}-h \nu_0}{h \nu_0}\right\rfloor
```

```
 n  Eₙ          ΔE
```

