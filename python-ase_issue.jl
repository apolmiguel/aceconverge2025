ENV["JULIA_PKG_USE_CLI_GIT"] = "true"

import Pkg
Pkg.activate(".")
# Pkg.add(Pkg.PackageSpec(;name="ACEpotentials", version="0.6.11"))
using ACEpotentials, LinearAlgebra, Plots, LaTeXStrings

# dataset and model 
rawdata, _, _ = ACEpotentials.example_dataset("Si_tiny")
datakeys = (energy_key = "dft_energy", force_key)

# dataset
data_file = "/leonardo_work/Sis25_degironc_0/apol/Tr124_dim.xyz"
data = read_extxyz(data_file);