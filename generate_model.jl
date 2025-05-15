ENV["JULIA_PKG_USE_CLI_GIT"] = "true"
cd("/leonardo_work/Sis25_degironc_0/apol/aceconverge2025/")
pwd()

import Pkg
Pkg.activate(".")

using ACEpotentials, DelimitedFiles

# Dimer exporter
function export_dimers_to_dat(dimers_dict; filename)
    for ((z1, z0), (rr, v01)) in dimers_dict
        open(filename, "w") do io
            for (r, v) in zip(rr, v01)
                println(io, "$r $v")
            end
        end
    end
end

# Load data
traindir = ARGS
traindir = "/leonardo_work/Sis25_degironc_0/apol/aceconverge2025/datasets/Tr124_dim.xyz"
train_data = read_extxyz(traindir);
valdir = "/leonardo_work/Sis25_degironc_0/apol/aceconverge2025/datasets/Val123_dim.xyz"
val_data = read_extxyz(valdir);



# # basis construction with purification 
r0 = 1.5
basis_vanilla = ACE1x.ace_basis(
    elements = [:C],
    order = 3,
    totaldegree = [45, 5, 4] , # issue with setting 
    rcut = 7.0,
    r0 = r0,
    pure = false
)
# basis_purified = ACE1x.ace_basis(
#     elements = [:C],
#     order = 3,
#     totaldegree = [45, 5, 4] , # issue with setting 
#     rcut = 7.0,
#     r0 = r0,
#     pure = true
# )


println("Length of basis is: ", length(basis_vanilla))

println("\nAssigning offset")
Vref = OneBody(:C => -245.44385736) # one-body energy
println("Vref: ", Vref)
println("\nAssigning weights")
weights = Dict("shaiducarbon" => Dict("E" => 50.0, "F" => 1.0)) # loss weights
println("Weights: ", weights)

# Basis precomputation 
datakeys = (energy_key = "energy", force_key = "force")
println("\nConstructing linear problem elements")
train_atoms = [ACEpotentials.AtomsData(t; weights=weights, v_ref=Vref, datakeys...) for t in train_data]
A, Y, W = ACEfit.assemble(train_atoms, basis_vanilla)
println("Creating prior")
P = smoothness_prior(basis_vanilla; p=2) # higher p, stronger regularization 

# Potential 1
println("\nAssigning solver with prior")
solver = ACEfit.LSQR(damp = 1e-4, atol = 1e-6, P=P)
println("Solving linear problem")
results = ACEfit.solve(solver, W .* A, W .* Y)
pot_1 = JuLIP.MLIPs.SumIP(Vref, JuLIP.MLIPs.combine(basis_vanilla, results["C"]))

# # Potential 2
# # solver = ACEfit.LSQR(damp = 1e-4, atol = 1e-6, P=P)
# # results = ACEfit.solve(solver, W .* A, W .* Y)

# Errors and validation 
println("\nValidating")
val_atoms = [ACEpotentials.AtomsData(t; weights=weights, v_ref=Vref, datakeys...) for t in val_data]
err1 = ACEpotentials.linear_errors(train_atoms, pot_1)
println("Saving errors to file datafiles/err1.dat")
writedlm("datafiles/err1.dat", err1)

# Dimers
D = ACEpotentials.dimers(pot_1, [:C,]; rr = range(0.529177, 7.0, length=200))
export_dimers_to_dat(D, filename="datafiles/dim1.dat")

