import Pkg
Pkg.activate(".")

using ACEpotentials, DelimitedFiles

## Functions ## 
function export_dimers_to_dat(dimers_dict; filename)
    for ((z1, z0), (rr, v01)) in dimers_dict
        open(filename, "w") do io
            for (r, v) in zip(rr, v01)
                println(io, "$r $v")
            end
        end
    end
end

## Data load ##
traindir = ARGS[1]
train_data = read_extxyz(traindir);
valdir = ARGS[2]
val_data = read_extxyz(valdir);
pureflag = ARGS[3]
elossweight = parse(Float64, ARGS[4]) 

orders = [2,3,4]
degrees = [[24,20],[24,20,16],[24,20,16,12]]
basis_tags = ["24.20","24.20.16","24.20.16.12"]

# orders = [2,3,3,3,4]
# degrees = [[46,16],[46,16,12],[46,20,14],[46,24,16],[46,20,14,10]]
# basis_tags  = ["46.16","46.16.12","46.20.14","46.24.16","46.20.14.10"]

## Assigning common solver properties from training ## 
println("\nAssigning offset")
Vref = OneBody(:C => -245.44385736) # one-body energy
println("Vref: ", Vref)
println("\nAssigning weights")
weights = Dict("shaiducarbon" => Dict("E" => elossweight, "F" => 1.0)) # loss weights
println("Weights: ", weights)
datakeys = (energy_key = "energy", force_key = "force")
train_atoms = [ACEpotentials.AtomsData(t; weights=weights, v_ref=Vref, datakeys...) for t in train_data]
val_atoms = [ACEpotentials.AtomsData(t; weights=weights, v_ref=Vref, datakeys...) for t in val_data]


if ARGS[3] == "purify"
    pureflag = true # canonical CE
    prefix = "acejulia/" * traindir[10:end-4] * "_purify/" # for saving files
else
    pureflag = false # self-interacting CE
    prefix = "acejulia/" * traindir[10:end-4] * "/" # for saving files
end

for (i, label) in enumerate(basis_tags)
    rundir = prefix * label * "/" * "ecost$(elossweight)/"
    potdir = rundir * "potential.json"
    errdir = rundir * "errors.dat"
    dimerdir = rundir * "dimercurve.dat"
    println("Loading potential from $(potdir).")
    pot = load_potential(potdir)
    err = ACEpotentials.linear_errors(train_atoms, pot)
    println("Writing errors to $(errdir).")
    writedlm(errdir, err)
    D = ACEpotentials.dimers(pot, [:C,]; rr = range(0.529177, 7.0, length=200))
    println("Writing dimer curves to $(dimerdir).")
    export_dimers_to_dat(D, filename=dimerdir)
end


