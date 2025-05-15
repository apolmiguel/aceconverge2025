import Pkg
Pkg.activate(".")
# Pkg.status()

using ACEpotentials, DelimitedFiles


## Variables ## 
# user-input
rcut = 7.0 # 7.0 if dimer, 5.0 if diamond or generic
pureflag = false # true if want to get canonical CE

# controlled
orders = [2,3]
degrees = [[40,10], [40,10,9]] # actual
# degrees = [[5,5], [5,5,5]] # testing 
r0 = 1.286958464 # minimum from dimer dataset 

# create bases 
basis_bin = Dict()
for (i,ord) in enumerate(orders)
    print("i = $i, ord = $ord, degrees = $(degrees[i])\n")
    basis = ACE1x.ace_basis(
        elements = [:C],
        order = ord,
        totaldegree = degrees[i],
        rcut = rcut,
        r0 = r0,
        pure = pureflag
    )
    basis_bin["b_order$(ord+1)"] = basis
    println("basis_bin[\"b_order$(ord+1)\"] = ", basis_bin["b_order$(ord+1)"])
end


## Data load ## 
# traindir = "/leonardo_work/Sis25_degironc_0/apol/aceconverge2025/datasets/Tr124_dim.xyz"
traindir = ARGS[1]
train_data = read_extxyz(traindir);
# valdir = "/leonardo_work/Sis25_degironc_0/apol/aceconverge2025/datasets/Val123_dim.xyz"
valdir = ARGS[2]
val_data = read_extxyz(valdir);
println("train_data: ", length(train_data), " val_data: ", length(val_data))



# println("Length of basis is: ", length(basis_vanilla))

# println("\nAssigning offset")
# Vref = OneBody(:C => -245.44385736) # one-body energy
# println("Vref: ", Vref)
# println("\nAssigning weights")
# weights = Dict("shaiducarbon" => Dict("E" => 50.0, "F" => 1.0)) # loss weights
# println("Weights: ", weights)

# # Basis precomputation 
# datakeys = (energy_key = "energy", force_key = "force")
# println("\nConstructing linear problem elements")
# train_atoms = [ACEpotentials.AtomsData(t; weights=weights, v_ref=Vref, datakeys...) for t in train_data]
# A, Y, W = ACEfit.assemble(train_atoms, basis_vanilla)
# println("Creating prior")
# P = smoothness_prior(basis_vanilla; p=2) # higher p, stronger regularization 

# # Potential 1
# println("\nAssigning solver with prior")
# solver = ACEfit.LSQR(damp = 1e-4, atol = 1e-6, P=P)
# println("Solving linear problem")
# results = ACEfit.solve(solver, W .* A, W .* Y)
# pot_1 = JuLIP.MLIPs.SumIP(Vref, JuLIP.MLIPs.combine(basis_vanilla, results["C"]))

# # # Potential 2
# # # solver = ACEfit.LSQR(damp = 1e-4, atol = 1e-6, P=P)
# # # results = ACEfit.solve(solver, W .* A, W .* Y)

# # Errors and validation 
# println("\nValidating")
# val_atoms = [ACEpotentials.AtomsData(t; weights=weights, v_ref=Vref, datakeys...) for t in val_data]
# err1 = ACEpotentials.linear_errors(train_atoms, pot_1)
# println("Saving errors to file datafiles/err1.dat")
# writedlm("datafiles/err1.dat", err1)

# # Dimers
# D = ACEpotentials.dimers(pot_1, [:C,]; rr = range(0.529177, 7.0, length=200))
# export_dimers_to_dat(D, filename="datafiles/dim1.dat")

