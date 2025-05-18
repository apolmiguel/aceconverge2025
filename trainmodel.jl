import Pkg
Pkg.activate(".")
using ACEpotentials, DelimitedFiles

## Explicit parsing ## 
traindir = ARGS[1] # training data load
println("Training dataset taken from $traindir.")
train_data = read_extxyz(traindir);
if ARGS[2] == "dim" # rcut flag
    rcut = 7.0 
else
    rcut = 5.0 # dia, gen cutoff
end
println("rcut = $rcut")
if ARGS[3] == "purify" # purification flag
    pureflag = true 
    prefix = "acejulia/" * traindir[10:end-4] * "_purify/"
else
    pureflag = false # self-interacting CE
    prefix = "acejulia/" * traindir[10:end-4] * "/" 
end
println("Are we using the purification procedure: $pureflag.")
elossweight= parse(Float64, ARGS[4]) # energy weight
println("(with Fcost = 1.0), Ecost = $elossweight.")
dampval = parse(Float64, ARGS[5]) # solver damping value
println("Damping value for LSQR solver = $dampval.")
solverflag = ARGS[6]
println("Solver flag: $solverflag.")


## Control parameters ##
orders = [2,2,3,4]
degrees = [[10,5], [40,10], [40,10,9], [40,10,9,8]]
basis_tags  = ["10.5","40.10","40.10.9", "40.10.9.8"]

# orders = [2,3,3,3,4]
# degrees = [[46,16],[46,16,12],[46,20,14],[46,24,16],[46,20,14,10]]
# basis_tags  = ["46.16","46.16.12","46.20.14","46.24.16","46.20.14.10"]
r0 = 1.286958464 # equilibrium length from dimer dataset


## Common solver properties ## 
println("\nAssigning offset.")
Vref = OneBody(:C => -245.44385736) # one-body energy
println("Vref: ", Vref)
println("Assigning weights.")
weights = Dict("shaiducarbon" => Dict("E" => elossweight, "F" => 1.0))
println("Weights: ", weights)
datakeys = (energy_key = "energy", force_key = "force")
train_atoms = [ACEpotentials.AtomsData(t; weights=weights, v_ref=Vref, datakeys...) for t in train_data]
# basis_bin = Dict()


for (i, label) in enumerate(basis_tags)
    # Basis creation
    println("\nCreating basis for order $(orders[i]), with per-correlation degrees $(degrees[i]).")
    println("r0 = $r0, rcut = $rcut.")
    basis = ACE1x.ace_basis(
        elements = [:C],
        order = orders[i],
        totaldegree = degrees[i],
        rcut = rcut,
        r0 = r0,
        pure = pureflag)
    println("Basis for $label created, with $(length(basis)) basis functions.")

    # Linear problem assembly
    potdir = prefix * label * "/" * "ecost$(elossweight)/" * "damp$(dampval)/"
    println("Creating directory: $potdir")
    mkpath(potdir)

    println("\nAssembling linear problem elements: A, Y, W for basis set: $label.")
    A, Y, W = ACEfit.assemble(train_atoms, basis)
    println("Creating prior.")
    P = smoothness_prior(basis; p=2)
    println("Creating solver.")
    if length(ARGS) >= 6 && solverflag == "BLR"
        solver = ACEfit.BLR()
    else
        solver = ACEfit.LSQR(damp = dampval, atol = 1e-6, P = P)
    end
    println("Solving linear problem.")
    results = ACEfit.solve(solver, W .* A, W .* Y)
    println("Creating potential.")
    pot = JuLIP.MLIPs.SumIP(Vref, JuLIP.MLIPs.combine(basis, results["C"]))
    if solverflag == "BLR"
        pot = ACEpotentials.BLRIP(pot, results["C"], results["Cerr"])
        println("Saving potential at $potdir.")
        save_potential(potdir * "potential_BLR.json", pot)
    else
        println("Saving potential at $potdir.")
        save_potential(potdir * "potential.json", pot)
    end
end
