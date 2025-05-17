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
# elossweight = parse(Float64, ARGS[4])
elossweight = 50.0
# dampval = parse(Float64, ARGS[5])
dampval = [1e-15, 1e-9, 5e-9, 1e-8]
## Control parameters ## 
orders = [2,3,3,3,4]
degrees = [[46,16],[46,16,12],[46,20,14],[46,24,16],[46,20,14,10]]
basis_tags  = ["46.16","46.16.12","46.20.14","46.24.16","46.20.14.10"]

if ARGS[3] == "purify"
    pureflag = true # canonical CE
    prefix = "acejulia/" * traindir[10:end-4] * "_purify/" # for saving files
else
    pureflag = false # self-interacting CE
    prefix = "acejulia/" * traindir[10:end-4] * "/" # for saving files
end

for (i, label) in enumerate(basis_tags)
    for (j, damp) in enumerate(dampval)
        rundir = prefix * label * "/" * "ecost$(elossweight)/" * "damp$(damp)/"
        potdir = rundir * "potential.json"
        errdir = rundir * "errors.dat"
        dimerdir = rundir * "dimercurve.dat"
        println("Loading potential from $(potdir).")
    end
end



# ## Assigning common solver properties from training ## 
# println("\nAssigning offset")
# Vref = OneBody(:C => -245.44385736) # one-body energy
# println("Vref: ", Vref)
# println("\nAssigning weights")
# weights = Dict("shaiducarbon" => Dict("E" => 50.0, "F" => 1.0)) # loss weights
# println("Weights: ", weights)
# datakeys = (energy_key = "energy", force_key = "force")
# train_atoms = [ACEpotentials.AtomsData(t; weights=weights, v_ref=Vref, datakeys...) for t in train_data]

# ## Validation ## 
# val_atoms = [ACEpotentials.AtomsData(t; weights=weights, v_ref=Vref, datakeys...) for t in val_data]

# # prefix = "acejulia/" * traindir[10:end-4]


# for j in 3:5
#     # Directories dependent on purification 
#     # if pureflag == "purify"
#     #     potdir = prefix * "border$(j)/potential_pure.json"
#     #     errdir = prefix * "border$(j)/errors_pure.dat"
#     #     dimerdir = prefix * "border$(j)/dimercurve_pure.dat"
#     # else
#     potdir = prefix * "border$(j)/potential.json"
#     errdir = prefix * "border$(j)/errors.dat"
#     dimerdir = prefix * "border$(j)/dimercurve.dat"
#     # end
    
#     # Load potential
#     println("Loading potential from $(potdir).")
#     pot = load_potential(potdir)

#     # Error analysis
#     err = ACEpotentials.linear_errors(train_atoms, pot)
#     println("Writing errors to $(errdir).")
#     writedlm(errdir, err)

#     # Dimer curves\
#     D = ACEpotentials.dimers(pot, [:C,]; rr = range(0.529177, 7.0, length=200))
#     println("Writing dimer curves to $(dimerdir).")
#     export_dimers_to_dat(D, filename=dimerdir)
# end