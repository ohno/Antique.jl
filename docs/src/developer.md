# [Developer Guide](@id developer-guide)

This is the guideline for adding new models. Adding a new model may take from a few days to a few weeks due to reference search, test implementation, and writing documentation.

1. First, please submit a new issue or comment [here](https://github.com/ohno/Antique.jl/issues). I will assign you to the issue. We need to find orthodox references (textbooks or papers, not Wikipedia) for the analytical solutions (eigenvalues and eigenfunctions) before the development. This will take more time than you think.
2. Fork [the repository](https://github.com/ohno/Antique.jl) on GitHub.
3. Clone the forked repository to your local machine by Git.
4. Please create 3 files:

| files | comments |
| --- | --- |
| `src/ModelName.jl` | Write the source codes and docstrings in this file. The most helpful examples are the harmonic oscillator for one-dimensional systems and the hydrogen atom for three-dimensional systems. We recommend that you copy these files. First create a structure `struct ModelName` with the same name as the model name (The best way is Find & Replace). Because the function names conflict, you must always give the struct `ModelName` as the first argument to V, E, ψ and other functions. Multi-dispatch avoids conflicts. We recommend using Revice.jl while coding. Run `include("./dev/revice.jl")` on the REPL or use dev.ipynb. |
| `test/ModelName.jl` | Write test code in this file. At a minimum, please check the normalization and the orthogonality of the eigenfunctions using QuadGK.jl. Please also do tests for the eigenvalues (for example, calculate the expectation values of the Hamiltonian (energy) using the eigenfunctions and check that these values match the eigenvalues). |
| `docs/src/ModelName.md` | Write documentation in this file. Include at least the definition of the Hamiltonian and the analytical solutions (eigenvalues and eigenfunctions). Call a docstring in the source code (`src/ModelName.jl`) . |

5. Please rewrite 5 files:

| files | comments |
| - | - |
| `src/Antique.jl` | Add the new model name `:ModelName` to the `models = [...]` array in this file. `:` is required at the beginning. |
| `docs/make.jl` | Add the new model into `pages=[...]` in this file. |
| `test/runtests.jl` | Change `for model in [...]` in this file. Please test all models before pull requests. |
| `README.md` | Add the new model to the list of supported models. |
| `docs/index.md` | Add the new model to the list of supported models. |

6. Execute `include("./dev/test.jl")` to run tests. It will take few minutes to complete.
7. Execute `include("./dev/docs.jl")` to compile documents. HTML files (docs/build/*.html) will be generated. Please check them with Chrome or any other web browsers.
8. Commit and Push the codes.
9. Submit a pull request on GitHub.

## Acknowledgment

Thanks to all contributors. This package was named by [@KB-satou](https://github.com/KB-satou) and [@ultimatile](https://github.com/ultimatile). [@MartinMikkelsen](https://github.com/MartinMikkelsen) contributed to writing docstrings. Special thanks to [@hyrodium](https://github.com/hyrodium) for his help with managing the documentation and advice on coding style. [@lhapp27](https://github.com/lhapp27) implemented 2 models, and [@ajarifi](https://github.com/ajarifi) implemented 3 models.
