import Pkg
Pkg.activate(".")
using ACEpotentials, DelimitedFiles

## Collect CLI arguments ##
basis_tag = ARGS[1]
println("Basis tag: $basis_tag")