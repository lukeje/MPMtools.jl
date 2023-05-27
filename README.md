# MPMtools.jl
Some useful scripts for [MPM](https://doi.org/10.1016/j.neuroimage.2019.01.029) MRI data

## Installation
This project is not currently listed in the `julia` package management system, and thus must be manually installed.
One way to do this is using an environment.

First clone the repository on the command line:
```bash
git clone https://github.com/lukeje/MPMtools.jl/
julia # open julia REPL
```
then in the `julia` REPL run
```julia
using Pkg
Pkg.activate("MPMtools.jl")
Pkg.instantiate()
```

The environment can then be used within `julia` by running
```julia
using Pkg
Pkg.activate("MPMtools.jl")
using MPMtools
```

## Requirements
An installation of a recent version of `julia`.

## Current status
This repository is a result of me wanting to try out the [`julia` language](https://docs.julialang.org/).
As I'm still new to the language, I make no guarantees that the scripts work as expected!