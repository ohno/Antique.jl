#### Normalization & Orthogonality of $\psi_n(x)$
```math
\int_{0}^{L} \psi_i^\ast(x) \psi_j(x) ~\mathrm{d}x = \delta_{ij}
```

```
 i |  j |     analytical |      numerical 
-- | -- | -------------- | -------------- 
```

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

```
  L |   m |   ℏ |  n |     analytical |      numerical 
--- | --- | --- | -- | -------------- | -------------- 
```

#### Expected Value of $x$

```math
\langle x \rangle_{n=1}
= \int_{0}^{L} \psi_1^\ast(x) \hat{x} \psi_1(x) ~\mathrm{d}x
= \frac{2(2a)^2}{\pi^3} \left( \frac{\pi^3}{6} - \frac{\pi}{4} \right)
```
for only $n=1$.

Reference:
- [LibreTexts PHYSICS, 6.4: Expectation Values, Observables, and Uncertainty](https://phys.libretexts.org/Bookshelves/Modern_Physics/Book%3A_Spiral_Modern_Physics_(D'Alessandris)/6%3A_The_Schrodinger_Equation/6.4%3A_Expectation_Values_Observables_and_Uncertainty)
```
  L |  n |     analytical |      numerical 
--- | -- | -------------- | -------------- 
```

#### Expected Value of $x^2$

```math
\langle x^2 \rangle_{n=1}
= \int_{0}^{L} \psi_1^\ast(x) \hat{x}^2 \psi_1(x) ~\mathrm{d}x
= \frac{2(2a)^2}{\pi^3} \left( \frac{\pi^3}{6} - \frac{\pi}{4} \right)
```

Reference:
- [LibreTexts PHYSICS, 6.4: Expectation Values, Observables, and Uncertainty](https://phys.libretexts.org/Bookshelves/Modern_Physics/Book%3A_Spiral_Modern_Physics_(D'Alessandris)/6%3A_The_Schrodinger_Equation/6.4%3A_Expectation_Values_Observables_and_Uncertainty)
```
  L |  n |     analytical |      numerical 
--- | -- | -------------- | -------------- 
```


#### Expected Value of $p$
```math
\langle p \rangle_{n=1}
= \int_{0}^{L} \psi_1^\ast(x) \hat{p} \psi_1(x) ~\mathrm{d}x
= 0
```

Reference:
- [LibreTexts PHYSICS, 6.4: Expectation Values, Observables, and Uncertainty](https://phys.libretexts.org/Bookshelves/Modern_Physics/Book%3A_Spiral_Modern_Physics_(D'Alessandris)/6%3A_The_Schrodinger_Equation/6.4%3A_Expectation_Values_Observables_and_Uncertainty)

---

```math
  \begin{aligned}
    \langle p \rangle_{n=1}
    &= \int_0^L \psi^\ast_n(x) \hat{p} \psi_n(x) ~\mathrm{d}x \\
    &= \int_0^L \psi^\ast_n(x) \left[ -i\hbar\frac{\mathrm{d}}{\mathrm{d} x} \right] \psi(x) ~\mathrm{d}x \\
    &\simeq \int_0^L \psi^\ast_n(x) \left[ -i\hbar \frac{\psi(x+\Delta x) - \psi(x-\Delta x)}{2\Delta x} \right] ~\mathrm{d}x.
  \end{aligned}
```

Where, the difference formula for the 2nd-order derivative:

```math
\begin{aligned}
  % 2\frac{\mathrm{d} \psi(x)}{\mathrm{d}x} \Delta x
  % + O\left(\Delta x^{3}\right)
  % &= 
  % \psi(x+\Delta x)
  % - \psi(x-\Delta x)
  % \\
  % 2\frac{\mathrm{d} \psi(x)}{\mathrm{d}x} \Delta x
  % &= 
  % \psi(x+\Delta x)
  % - \psi(x-\Delta x)
  % - O\left(\Delta x^{3}\right)
  % \\
  % \frac{\mathrm{d} \psi(x)}{\mathrm{d}x}
  % &= 
  % \frac{\psi(x+\Delta x)- \psi(x-\Delta x)}{2\Delta x}
  % - \frac{O\left(\Delta x^{3}\right)}{2\Delta x}
  % \\
  \frac{\mathrm{d} \psi(x)}{\mathrm{d}x}
  &= 
  \frac{\psi(x+\Delta x)- \psi(x-\Delta x)}{2\Delta x}
  + O\left(\Delta x^{2}\right),
\end{aligned}
```

are given by the sum of 2 Taylor series:

```math
\begin{aligned}
  \psi(x+\Delta x)
  &=
  \psi(x)
  + \frac{\mathrm{d} \psi(x)}{\mathrm{d}x} \Delta x
  + \frac{1}{2!} \frac{\mathrm{d}^{2} \psi(x)}{\mathrm{d}x^{2}} \Delta x^{2}
  + O\left(\Delta x^{3}\right),
  \\
  \psi(x-\Delta x)
  &=
  \psi(x)
  - \frac{\mathrm{d} \psi(x)}{\mathrm{d}x} \Delta x
  + \frac{1}{2!} \frac{\mathrm{d}^{2} \psi(x)}{\mathrm{d}x^{2}} \Delta x^{2}
  + O\left(\Delta x^{3}\right).
\end{aligned}
```

```
  L |  n |     analytical |      numerical 
--- | -- | -------------- | -------------- 
```

#### Expected Value of $p^2$

```math
\langle p^2 \rangle
= \int_{0}^{L} \psi_1^\ast(x) \hat{p}^2 \psi_1(x) ~\mathrm{d}x
= \frac{\pi^2\hbar^2}{L^2}
```

Reference:
- [LibreTexts PHYSICS, 6.4: Expectation Values, Observables, and Uncertainty](https://phys.libretexts.org/Bookshelves/Modern_Physics/Book%3A_Spiral_Modern_Physics_(D'Alessandris)/6%3A_The_Schrodinger_Equation/6.4%3A_Expectation_Values_Observables_and_Uncertainty)

---

```math
  \begin{aligned}
    \langle p^2 \rangle
    &= \int_0^L \psi^\ast_n(x) \hat{p} \psi_n(x) ~\mathrm{d}x \\
    &= \int_0^L \psi^\ast_n(x) \left[ -\hbar^2\frac{\mathrm{d}^2}{{\mathrm{d}x}^2} \right] \psi(x) ~\mathrm{d}x \\
    &\simeq \int_0^L \psi^\ast_n(x) \left[ -\hbar^2 \frac{\psi(x+\Delta x) - 2\psi(x) + \psi(x-\Delta x)}{\Delta x^{2}} \right] ~\mathrm{d}x.
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
  L |  n |     analytical |      numerical 
--- | -- | -------------- | -------------- 
```

