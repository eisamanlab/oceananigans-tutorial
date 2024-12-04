# Oceananigans.jl Tutorial

The purpose of this repository to get comfortable with [Oceananigans.jl](https://clima.github.io/OceananigansDocumentation/stable/).

# Getting Started

1. **Download Julia**

Download Julia from https://julialang.org/downloads/

2. **Activate and instantitae the environment**

Make sure you are in the directory containing the `Project.toml` and `Manifest.toml` files. Then, in the Julia REPL run the following
```julia
Using Pkg; # loads the Julia package manager
Pkg.activate("."); # activates the environment defined by Project.toml
Pkg.instantiate(); # installs dependencies defined by Manifest.toml
```
3. **Create a Julia Jupyter kernel**

From the Julia REPL, run the following:

```julia
using IJulia;
installkernel("julia"); # creates a julia kernel for Jupyter
```

4. **Open jupyterlab**

From the Julia REPL, run the following

```julia
using IJulia; # skip this if IJulia is already loaded
jupyterlab() # starts a jupyterlab server in this directory
```

5. **Start using Oceananigans!**

Now you can create a notebook and start experimenting with `Oceananigans.jl`. I suggest trying to run some of the examples in the [documentation](https://clima.github.io/OceananigansDocumentation/stable/).


 # Additional resources
  - [Julia Academy](https://juliaacademy.com/) is a great resource is you are just learning Julia
  - The Julia [Getting Started](https://docs.julialang.org/en/v1/manual/getting-started/) guide is another great resource
  