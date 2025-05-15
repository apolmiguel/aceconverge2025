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
potdir = ARGS[3]
potential = load_potential(potdir, verbose = true)

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

for j in 3:4
    # Error analysis
    err = ACEpotentials.linear_errors(train_atoms, potential)
    errdir = "datafiles/Tr124_dim_errors_b_order$(j).dat"
    println("Writing errors to $(errdir).")
    writedlm(errdir, err)

    # Dimer curves
    D = ACEpotentials.dimers(potential, [:C,]; rr = range(0.529177, 7.0, length=200))
    dimerdir = "datafiles/Tr124_dim_dimcurve_b_order$(j).dat"
    println("Writing dimer curves to $(dimerdir).")
    export_dimers_to_dat(D, filename=dimerdir)