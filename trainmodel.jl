import Pkg
Pkg.activate(".")

using ACEpotentials, DelimitedFiles


## Data load ## 
traindir = ARGS[1]
train_data = read_extxyz(traindir);
valdir = ARGS[2]
val_data = read_extxyz(valdir);
println("Training set has ", length(train_data), " entries.\nValidation set has ", length(val_data), " entries.")


## Variables ##
if ARGS[3] == "dim"
    rcut = 7.0 # cutoff for dimer dataset
else
    rcut = 5.0 # cutoff for diamond or generic 
end
println("rcut = $rcut")

if length(ARGS) >= 4 && ARGS[4] == "purify"
    pureflag = true # canonical CE
    prefix = "acejulia/" * traindir[10:end-4] * "_purify/" # for saving files

else
    pureflag = false # self-interacting CE
    prefix = "acejulia/" * traindir[10:end-4] * "/" # for saving files
end

# prefix = "acejulia/" * traindir[10:end-4] # for saving files 
println("pureflag = $pureflag")
orders = [2,3,3,3,4]
degrees = [[46,16],[46,16,12],[46,20,14],[46,24,16],[46,20,14,10]]
r0 = 1.286958464 # minimum from dimer dataset
println("r0 = $r0")


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

for j in 3:5 # body/correlation order = order + 1
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
    # if pureflag
    #     potdir = prefix * "/border$(j)/potential_pure.json"
    # else
    #     potdir = prefix * "/border$(j)/potential.json"
    # end
    potdir = prefix * "border$(j)/potential.json"
    println("Saving potential to file $potdir")
    save_potential(potdir, pot)
end
