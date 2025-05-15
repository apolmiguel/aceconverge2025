import Pkg
Pkg.activate(".")
# Pkg.status()

using ACEpotentials, DelimitedFiles


## Data load ## 
# traindir = "/leonardo_work/Sis25_degironc_0/apol/aceconverge2025/datasets/Tr124_dim.xyz"
traindir = ARGS[1]
train_data = read_extxyz(traindir);
# valdir = "/leonardo_work/Sis25_degironc_0/apol/aceconverge2025/datasets/Val123_dim.xyz"
valdir = ARGS[2]
val_data = read_extxyz(valdir);
println("Training set has ", length(train_data), " entries.\nValidation set has ", length(val_data), " entries.")


## Variables ##
if length(ARGS) >= 3 && ARGS[3] == "purify"
    pureflag = true # canonical CE
else
    pureflag = false # self-interacting CE
end
prefix = traindir[10:end-4]

orders = [2,3]
degrees = [[40,10], [40,10,9]] # actual
# degrees = [[5,5], [5,5,5]] # testing 
r0 = 1.286958464 # minimum from dimer dataset
rcut = 7.0 # for dimer 


## Basis creation ## 
println("\nCreating basis")
basis_bin = Dict()
for (i,ord) in enumerate(orders)
    println("i = $i, ord = $ord, degrees = $(degrees[i])\n")
    basis = ACE1x.ace_basis(
        elements = [:C],
        order = ord,
        totaldegree = degrees[i],
        rcut = rcut,
        r0 = r0,
        pure = pureflag
    )
    basis_bin["border$(ord+1)"] = basis
    println("basis_bin[\"border$(ord+1)\"] created.\n")
end

## Assigning common solver properties ## 
println("\nAssigning offset")
Vref = OneBody(:C => -245.44385736) # one-body energy
println("Vref: ", Vref)
println("\nAssigning weights")
weights = Dict("shaiducarbon" => Dict("E" => 50.0, "F" => 1.0)) # loss weights
println("Weights: ", weights)
datakeys = (energy_key = "energy", force_key = "force")
println("\nConstructing linear problem elements")
train_atoms = [ACEpotentials.AtomsData(t; weights=weights, v_ref=Vref, datakeys...) for t in train_data]

for j in 3:4 # body/correlation order = order + 1
    println("\nAssembling linear problem: A, Y, W for basis_bin[\"border$j\"]")
    A, Y, W = ACEfit.assemble(train_atoms, basis_bin["border$j"])
    println("Creating prior")
    P = smoothness_prior(basis_bin["border$j"]; p=2) # higher p, stronger regularization 
    println("Creating solver")
    solver = ACEfit.BLR()
    println("Solving linear problem")
    results = ACEfit.solve(solver, W .* A, W .* Y)
    pot = JuLIP.MLIPs.SumIP(Vref, JuLIP.MLIPs.combine(basis_bin["border$j"], results["C"]))

    # Saving potential
    potdir = "acejulia/" * prefix * "/border$(j)/potential.json"
    println("Saving potential to file $potdir")
    save_potential(potdir, pot)
end

    # println("Validating")
    # val_atoms = [ACEpotentials.AtomsData(t; weights=weights, v_ref=Vref, datakeys...) for t in val_data]
    # err = ACEpotentials.linear_errors(train_atoms, pot)
    # println("Saving errors to file datafiles/Tr124_dim_errors_border$j.dat")
    # writedlm("datafiles/Tr124_dim_errors_border$j.dat", err)

    # # println("Dimer curves")
    # D = ACEpotentials.dimers(pot, [:C,]; rr = range(0.529177, 7.0, length=200))
    # export_dimers_to_dat(D, filename="datafiles/Tr124_dim__dimcurve_border$j.dat")
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

