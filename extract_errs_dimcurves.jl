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


## Assigning common solver properties from training ## 
println("\nAssigning offset")
Vref = OneBody(:C => -245.44385736) # one-body energy
println("Vref: ", Vref)
println("\nAssigning weights")
weights = Dict("shaiducarbon" => Dict("E" => 50.0, "F" => 1.0)) # loss weights
println("Weights: ", weights)
datakeys = (energy_key = "energy", force_key = "force")
train_atoms = [ACEpotentials.AtomsData(t; weights=weights, v_ref=Vref, datakeys...) for t in train_data]

## Validation ## 
val_atoms = [ACEpotentials.AtomsData(t; weights=weights, v_ref=Vref, datakeys...) for t in val_data]

prefix = "acejulia/" * traindir[10:end-4]
for j in 3:4
    # Load potential
    potdir = prefix * "/border$(j)/potential_pure.json"
    println("Loading potential from $(potdir).")
    pot = load_potential(potdir)

    # Error analysis
    errdir = prefix * "/border$(j)/errors_pure.dat"
    err = ACEpotentials.linear_errors(train_atoms, pot)
    println("Writing errors to $(errdir).")
    writedlm(errdir, err)

    # Dimer curves
    dimerdir = prefix * "/border$(j)/dimercurve_pure.dat"
    D = ACEpotentials.dimers(pot, [:C,]; rr = range(0.529177, 7.0, length=200))
    println("Writing dimer curves to $(dimerdir).")
    export_dimers_to_dat(D, filename=dimerdir)
end