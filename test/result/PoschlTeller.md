#### Associated Legendre Polynomials $P_n^m(x)$

```math
  \begin{aligned}
    P_n^m(x)
    &= \left( 1-x^2 \right)^{m/2} \frac{\mathrm{d}^m}{\mathrm{d}x^m} P_n(x) \\
    &= \left( 1-x^2 \right)^{m/2} \frac{\mathrm{d}^m}{\mathrm{d}x^m} \frac{1}{2^n n!} \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left[ \left( x^2-1 \right)^n \right] \\
    &= \frac{1}{2^n} (1-x^2)^{m/2} \sum_{j=0}^{\left\lfloor\frac{n-m}{2}\right\rfloor} (-1)^j \frac{(2n-2j)!}{j! (n-j)! (n-2j-m)!} x^{(n-2j-m)}.
  \end{aligned}
```

``n=0, m=0:`` ✔
```math
\begin{aligned}
  P_{0}^{0}(x)
    = 1
  &= 1 \\
  &= 1
\end{aligned}
```

``n=1, m=0:`` ✔
```math
\begin{aligned}
  P_{1}^{0}(x)
    = \frac{1}{2} \frac{\mathrm{d}}{\mathrm{d}x} \left( -1 + x^{2} \right)
  &= x \\
  &= x
\end{aligned}
```

