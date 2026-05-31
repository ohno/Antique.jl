#### Associated Legendre Polynomials $P_n^m(x)$

```math
  \begin{aligned}
    P_n^m(x)
    &= \left( 1-x^2 \right)^{m/2} \frac{\mathrm{d}^m}{\mathrm{d}x^m} P_n(x) \\
    &= \left( 1-x^2 \right)^{m/2} \frac{\mathrm{d}^m}{\mathrm{d}x^m} \frac{1}{2^n n!} \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left[ \left( x^2-1 \right)^n \right] \\
    &= \frac{1}{2^n} (1-x^2)^{m/2} \sum_{j=0}^{\left\lfloor\frac{n-m}{2}\right\rfloor} (-1)^j \frac{(2n-2j)!}{j! (n-j)! (n-2j-m)!} x^{(n-2j-m)}.
  \end{aligned}
```

#### Normalization & Orthogonality of $P_n^m(x)$

```math
\int_{-1}^{1} P_i^m(x) P_j^m(x) \mathrm{d}x = \frac{2(j+m)!}{(2j+1)(j-m)!} \delta_{ij}
```

```
 m |  i |  j |     analytical |      numerical 
-- | -- | -- | -------------- | -------------- 
```

#### Normalization & Orthogonality of $Y_{lm}(\theta,\varphi)$

```math
\int_0^{2\pi}
\int_0^\pi
Y_{lm}(\theta,\varphi)^* Y_{l'm'}(\theta,\varphi) \sin(\theta)
~\mathrm{d}\theta \mathrm{d}\varphi
= \delta_{ll'} \delta_{mm'}
```

```
l₁ | l₂ | m₁ | m₂ |     analytical |      numerical 
-- | -- | -- | -- | -------------- | -------------- 
```

#### Associated Laguerre Polynomials $L_n^{k}(x)$

```math
  \begin{aligned}
  L_n^{k}(x)
    &= \frac{\mathrm{d}^k}{\mathrm{d}x^k} L_n(x) \\
    &= \frac{\mathrm{d}^k}{\mathrm{d}x^k} \frac{1}{n!} \mathrm{e}^x \frac{\mathrm{d}^n}{\mathrm{d}x ^n} \left( \mathrm{e}^{-x} x^n \right) \\
    &= \sum_{m=0}^{n-k} (-1)^{m+k} \frac{n!}{m!(m+k)!(n-m-k)!} x^m \\
    &= (-1)^k L_{n-k}^{(k)}(x)
  \end{aligned}
```

#### Normalization & Orthogonality of $L_n^{k}(x)$

```math
\int_{0}^{\infty} \mathrm{e}^{-x} x^k L_i^k(x) L_j^k(x) \mathrm{d}x = \frac{i!}{(i-k)!} \delta_{ij}
```

Replace $n+k$ with $n$ for [the definition of Wolfram MathWorld](https://mathworld.wolfram.com/AssociatedLaguerrePolynomial.html).
```
 i |  j |  k |     analytical |      numerical 
-- | -- | -- | -------------- | -------------- 
```

#### Normalization of $R_{nl}(r)$

```math
\int |R_{nl}(r)|^2 r^2 \mathrm{d}r = 1
```

```
 n |  l |     analytical |      numerical 
-- | -- | -------------- | -------------- 
```

#### Expected Value of $r$

```math
\langle r \rangle
= \int r |R_{n_1 l_1}(r)|^2 r^2 \mathrm{d}r
= \frac{a_\mu}{2Z} \left[ 3n^2 - l(l+1) \right] \\
a_\mu = a_0 \frac{m_\mathrm{e}}{\mu} \\
\frac{1}{\mu} = \frac{1}{m_\mathrm{e}} + \frac{1}{m_\mathrm{p}}
```

Reference:
- [高柳和夫『朝倉物理学大系 11 原子分子物理学』(2000, 朝倉書店) pp.11-22](https://www.asakura.co.jp/detail.php?book_code=13681)
- [Quan­tum Me­chan­ics for En­gi­neers by Leon van Dom­me­len](https://web1.eng.famu.fsu.edu/~dommelen/quantum/style_a/nt_rsexp.html)

```
 n |  l |     analytical |      numerical 
-- | -- | -------------- | -------------- 
```

#### Expected Value of $r^2$

```math
\langle r^2 \rangle
= \int r^2 |R_{n_1 l_1}(r)|^2 r^2 \mathrm{d}r
= \frac{a_\mu^2}{2Z^2} n^2 \left[ 5n^2 + 1 - 3l(l+1) \right] \\
a_\mu = a_0 \frac{m_\mathrm{e}}{\mu} \\
\frac{1}{\mu} = \frac{1}{m_\mathrm{e}} + \frac{1}{m_\mathrm{p}}
```

Reference:
- [高柳和夫『朝倉物理学大系 11 原子分子物理学』(2000, 朝倉書店) pp.11-22](https://www.asakura.co.jp/detail.php?book_code=13681)
- [Quan­tum Me­chan­ics for En­gi­neers by Leon van Dom­me­len](https://web1.eng.famu.fsu.edu/~dommelen/quantum/style_a/nt_rsexp.html)
```
 n |  l |     analytical |      numerical 
-- | -- | -------------- | -------------- 
```

#### Virial Theorem

The virial theorem $2\langle T \rangle + \langle V \rangle = 0$ and the definition of Hamiltonian $\langle H \rangle = \langle T \rangle + \langle V \rangle$ derive $\langle H \rangle = \frac{1}{2} \langle V \rangle$ and $\langle H \rangle = -\langle T \rangle$.

```math
\frac{1}{2} \int \psi_n^\ast(x) V(x) \psi_n(x) \mathrm{d}x = E_n
```

```
Antique.HydrogenAtom(Z=1, m_e=1.0, a_0=1.0, E_h=1.0, hbar=1.0)
 n |          analytical |           numerical 
-- | ------------------- | ------------------- 
```

#### Normalization & Orthogonality of $\psi_n(r,\theta,\varphi)$

```math
\int \psi_i^\ast(r,\theta,\varphi) \psi_j(r,\theta,\varphi) r^2 \sin(\theta) \mathrm{d}r \mathrm{d}\theta \mathrm{d}\varphi = \delta_{ij}
```

```
n₁ | n₂ | l₁ | l₂ | m₁ | m₂ |     analytical |      numerical 
-- | -- | -- | -- | -- | -- | -------------- | -------------- 
