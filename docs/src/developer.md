# [Developer Guide](@id developer-guide)

If you are planning significant changes, please open an [issue](https://github.com/ohno/Antique.jl/issues) first. The [ColPrac](https://github.com/SciML/ColPrac) guidelines are recommended. For Julia package development basics, see:
- [How to develop a Julia package](https://julialang.org/contribute/developing_package/)
- [Pkg: Creating packages](https://pkgdocs.julialang.org/v1/creating-packages/)

## Local Setup

This procedure is required only once. Install Git and Julia on your local machine before starting.

1. Fork [the repository](https://github.com/ohno/Antique.jl) on GitHub.
2. Clone the forked repository. Replace `xxxxxx` with your GitHub username.
   ```sh
   git clone https://github.com/xxxxxx/Antique.jl.git
   cd Antique.jl
   ```
3. Install development tools: [Revise.jl](https://github.com/timholy/Revise.jl) and [Runic.jl](https://github.com/fredrikekre/Runic.jl).
   ```sh
   julia --startup-file=no -e 'import Pkg; Pkg.add("Revise")'
   julia --project=@runic --startup-file=no -e 'using Pkg; Pkg.add("Runic")'
   ```

## Development Flow

This is the typical workflow for making changes.

1. Create a branch for your changes. Replace `xxx` with the issue number (e.g. `issue/123`).
   ```sh
   git switch -c issue/xxx
   ```
2. Start an interactive session with [Revise.jl](https://github.com/timholy/Revise.jl).
   ```sh
   julia --startup-file=no -i -e 'using Revise; import Pkg; Pkg.activate("."); using Antique'
   ```
3. Change the source code. When making new functions or updating docstrings, refer to [Documenter: Adding docstrings](https://documenter.juliadocs.org/stable/man/guide/#Adding-Some-Docstrings).
4. If you need a new dependency, use `julia --project=. --startup-file=no -e 'import Pkg; Pkg.add("SomePackage"); Pkg.resolve(); Pkg.instantiate()'`. Replace `SomePackage` with the actual package name.
5. Format the source code with [Runic.jl](https://github.com/fredrikekre/Runic.jl).
   ```sh
   julia --project=@runic --startup-file=no -e 'using Runic; exit(Runic.main(ARGS))' -- --inplace .
   ```
6. Run the tests. It will take a few minutes.
   ```sh
   julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'
   ```
7. Build the documentation locally. HTML files (docs/build/*.html) will be generated. Please check them with Chrome or any other web browsers.
   ```sh
   julia --project=docs --startup-file=no -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate();'
   julia --project=docs --startup-file=no -e 'include("docs/make.jl")'
   ```
8. Commit and Push the codes (after steps 5–7 succeed).
   ```sh
   git add "path/to/changed/file"
   git commit -m "commit message"
   git push origin issue/xxx
   ```
9. Submit a pull request on GitHub.

## Adding New Models

This is the guideline for adding new models. Adding a new model may take from a few days to a few weeks due to reference search, test implementation, and writing documentation.

1. First, please submit a new issue or comment [here](https://github.com/ohno/Antique.jl/issues). I will assign you to the issue. We need to find orthodox references (textbooks or papers, not Wikipedia) for the analytical solutions (eigenvalues and eigenfunctions) before the development. This will take more time than you think.
2. Finish [Local Setup](@ref) and follow [Development Flow](@ref).
3. Please create 3 files:
    1. `src/ModelName.jl` — Write the source code and docstrings in this file. The most helpful examples are the harmonic oscillator for one-dimensional systems and the hydrogen atom for three-dimensional systems. We recommend that you copy these files. First create a structure `struct ModelName` with the same name as the model name (The best way is Find & Replace). Because the function names conflict, you must always give the struct `ModelName` as the first argument to potential, energy, wavefunction and other functions. Multi-dispatch avoids conflicts. We recommend using Revise.jl while coding ([Development Flow](@ref)).
    2. `test/ModelName.jl` — Write test code in this file. At a minimum, please check the normalization and the orthogonality of the eigenfunctions using QuadGK.jl. Please also do tests for the eigenvalues (for example, calculate the expectation values of the Hamiltonian (energy) using the eigenfunctions and check that these values match the eigenvalues).
    3. `docs/src/ModelName.md` — Write documentation in this file. Include at least the definition of the Hamiltonian and the analytical solutions (eigenvalues and eigenfunctions). Call a docstring in the source code (`src/ModelName.jl`).
4. Please rewrite 3 files:
    1. `src/Antique.jl` — Add `export ModelName`, add `:ModelName` to `models = [...]`, add `include("./ModelName.jl")`, and add `using .ModelNames: ModelName`.
    2. `docs/make.jl` — Add the new model into `pages=[...]` in this file.
    3. `docs/src/index.md` — Add the new model to the list of supported models.
5. As described in [Development Flow](@ref), please complete the following steps: run the tests, build and review the documentation, commit, push, and submit a pull request.

## Versioning and Registering (for Maintainers)

This project follows [Semantic Versioning](https://semver.org/). When bumping the version, update the version number in:

- [Project.toml](https://github.com/ohno/Antique.jl/blob/main/Project.toml#L4)
- [index.md](https://github.com/ohno/Antique.jl/blob/main/docs/src/index.md)

To register this package in the [General](https://github.com/JuliaRegistries/General) registry, install [Registrator](https://github.com/JuliaRegistries/Registrator.jl?tab=readme-ov-file#install-registrator) and use it via the [GitHub App](https://github.com/JuliaRegistries/Registrator.jl?tab=readme-ov-file#via-the-github-app).