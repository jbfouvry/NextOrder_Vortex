#=
 julia --threads 8 run/Run.jl --eps 0.01 --q 0.1 \
 --sigma_LOR3S 1.0 --J_LOR3S 0.5 --J_LOC 1.0 --sigma_LOC 0.2 --kmax 10
=#
############################################################
include("../src/Main.jl")
############################################################
# Table of the J for which we compute a flux
# We take some margin away from the border
# of the considered integration domain
# ATTENTION, we cast the array to type
const TAB_J = [J_LOC + ((i - 1 // 2) * SIGMA_LOC / 1000) for i=-999:1000]
const NB_J = length(TAB_J) # Total number of J for which a flux is computed
#####
const TAB_FLUX = zeros(Float64,NB_J) # Container for the flux calculation
#####
# Arrays of Float64 used for the dump
const TAB_J_DUMP = [Float64(TAB_J[iJ]) for iJ = 1:NB_J]
const TAB_FLUX_DUMP = zeros(Float64,NB_J)
############################################################
# Computing all the fluxes
############################################################
function calc!()
    #####
    Threads.@threads for iJ = 1:NB_J
        #####
        J_loc = TAB_J[iJ] # Current action
        #####
        TAB_FLUX[iJ] = FluxRESCALED(KMAX,J_loc)
        #####
    end
end
############################################################
# Performing the calculation
############################################################
function run!()
    @time calc!() # Computing the flux
    #####
    println("Flux [J=",TAB_J[1000],"] = ",TAB_FLUX[1000])
end
############################################################
run!() # Running

