---
title: 'Antique.jl: A Julia package on analytical solutions of quantum mechanical equations'
tags:
  - Julia
  - physics
  - chemistry
  - quantum mechanics
authors:
  - name: Shuhei Ohno
    orcid: 0009-0001-5222-9726
    equal-contrib: false
    affiliation: "1, 2"
  - name: Ahmad Jafar Arifi
    orcid: 0000-0002-9530-8993
    equal-contrib: false
    affiliation: "1, 3"
  - name: Lucas Happ
    orcid: 0000-0002-3893-279X
    equal-contrib: false
    affiliation: "1"
affiliations:
 - name: Few-body Systems in Physics Laboratory, RIKEN Nishina Center, Wako 351-0198, Japan
   index: 1
 - name: Graduate School of Nanobioscience, Yokohama City University, Yokohama 236-0027, Japan
   index: 2
 - name: Research Center for Nuclear Physics, Osaka University, Ibaraki 567-0047, Japan
   index: 3
date: 25 April 2025
bibliography: paper.bib
---

# Summary

Antique.jl provides implementations of analytical solutions to solvable quantum mechanical models. Analytical solutions are the most reliable benchmarks for software testing in the development of numerical methods. In addition to testing numerical methods, this package is useful for teaching quantum mechanics. Each solution has been tested using both a computer algebra system and numerical integration thanks to the multiple dispatch in the Julia programming language.

# Statement of need

In modern physics and chemistry, numerical methods for quantum mechanical equations have been developed to calculate the properties of hadrons, nuclei, atoms, molecules, and other systems. The quantum mechanical models with analytical solutions, such as harmonic oscillator and hydrogen atom, are often used as benchmarks for software testing in the development of numerical methods [@Hiyama2018]. There are many articles and books on analytical solutions, but few software implementations. To our knowledge, Antique.jl is the first package on analytical solutions of solvable quantum mechanical models. The implementation for each model provides functions for the potential, the eigenvalues, the eigenfunctions and some other properties. Ten models are currently supported and more models will be added in the future.

Antique.jl was developed for researchers, lecturers, students, and any person who is interested in quantum mechanics. Since Julia is designed from the ground up for numerical and scientific computing [@Bezanson2017], many packages for quantum mechanical calculations will be developed using this package in the future. In addition to benchmarking numerical methods, possible applications range from prototyping quark models in hadron physics to developing basis sets in quantum chemistry and teaching materials for quantum mechanics. The documentation includes graphs suitable for teaching materials. We hope that textbooks for quantum mechanics will be written using this package.

# State of the field

To our knowledge, there is no dedicated package on analytical solutions, only a few implementations are available [@Velieva2020; @RodrguezGmez2020]. However, there are several packages on the special polynomials contained in the analytical solutions, such as associated Laguerre polynomials and associated Legendre polynomials which appear in the wave functions of hydrogen-like atoms.

The special polynomials are usually defined by Rodrigues' formula, but their closed-form is better for a direct implementation. The correspondence between a Rodrigues' formula and a closed-form is not usually provided in textbooks. Therefore, a closed-form or an implementation corresponding to the Rodrigues' formula is needed. However, the definitions vary from textbook to textbook [@Greiner2001; @Griffiths2018], and package to package [@SciPy2020; @Cpp]. The time and effort to survey, implement, and test the special polynomials is a barrier to implementing analytical solutions.

The Julia language removes this barrier by the multiple dispatch. The multiple dispatch makes it easier to test polynomials using both computer algebra system [@Gowda2022] and numerical integration [@quadgk]. The match of the closed-form and the expansion of the Rodrigues' formula is checked using a computer algebra system. Additionally, the orthonormality has been tested using numerical integration. Finally, the orthonormality of the wave functions is tested using numerical integration, and the expectation value of the Hamiltonian is calculated from the eigenfunction to confirm its correspondence with the eigenvalues.

# Example usage

The quantum mechanical models with analytical solutions are useful for software testing in the development of numerical methods. Figure 1 shows an example of a variational calculation for the hydrogen atom [@Thijssen2007] and a comparison with the exact solution provided by Antique.jl.

![\label{fig:usage}The radial density of four s-wave states of the hydrogen atom calculated by the Rayleigh-Ritz method, implemented in TwoBody.jl [@TwoBody]. We employed Gaussian basis functions $\phi_n(r) = \exp(-\nu_n r^2)$ whose exponents were determined by the geometric progression defined in the previous study [@Hiyama2018]. The numerical solution (solid blue line) is compared to the analytical solution (dashed black line) using Antique.jl.](./figure.pdf)

# Acknowledgement

We acknowledge all contributors and sponsors. Special thanks to Yuto Horikawa for his help with managing the documentation and advice on the coding style. S. O. was supported by the RIKEN Junior Research Associate Program. A. J. A. and L. H. were supported by the RIKEN Special Postdoctoral Researcher Program.

# References
