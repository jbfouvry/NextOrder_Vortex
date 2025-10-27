#=
julia --threads 8 run/Run_Flux.jl --type Float64 --eps 0.01 \
--q 0.1 --sigma_LOR3 1.0 --J_LOC 1.0 --sigma_LOC 0.2 --nb_int 1000 --kmax 10
=#
############################################################
include("../src/Main.jl")
############################################################
# Table of the J for which we compute a flux
# We take some margin away from the border
# of the considered integration domain
# ATTENTION, we cast the array to type
const TAB_J = [TYPE(J_LOC + ((i - 1 // 2) * SIGMA_LOC / 100)) for i=-99:100]
const NB_J = length(TAB_J) # Total number of J for which a flux is computed
#####
const TAB_FLUX = zeros(TYPE,NB_J) # Container for the flux calculation
#####
# Arrays of Float64 used for the dump
const TAB_J_DUMP = [Float64(TAB_J[iJ]) for iJ = 1:NB_J]
const TAB_FLUX_DUMP = zeros(Float64,NB_J)
############################################################
# Performs the calculations over all the fundamental resonances
# that satisfy max(k,kp) = KMAX
# @IMPROVE -- The test to spot the fundamental resonances is very lazy
############################################################
function calc!()
    for k=1:KMAX, kp=1:KMAX # Loop over the fundamental resonances
        if (max(k,kp) == KMAX)
            calc_fund!(k,kp) # Computing for the current fundamental resonance
        end
    end
    #####
    # Casting the result to Float64
    for iJ = 1:NB_J
        TAB_FLUX_DUMP[iJ] = Float64(TAB_FLUX[iJ])
    end
end
############################################################
# Calculation over a given fundamental resonance
############################################################
function calc_fund!(k::Int64,kp::Int64)
    Threads.@threads for iJ=1:NB_J # Loop over the considered actions
        J_loc = TAB_J[iJ] # Current action
        #####
        TAB_FLUX[iJ] += FluxRESCALED(k,kp,J_loc) # Adding the contribution
    end
end
############################################################
# Performing the calculation
############################################################
function run!()
    @time calc!() # Computing the flux
    #####
    println("Flux [J=",TAB_J[20],"] = ",TAB_FLUX[20])
end
############################################################
run!() # Running
