#=
julia run/Run.jl --Npart 100 --q 0.0006 --eps 0.01 \
 --background LOR3S --recenter false --sigma_LOR3S 1.0 --J_LOR3S 0.5 \
 --J_LOC 1.0 --sigma_LOC 0.2 \
 --scheme RK9 --dt 2.0 --Nsteps 100 --Ndumps 1 --seed 1
=#
##################################################
include("../src/Main.jl")
##################################################
init!() # Initialising the particles
##################################################
# Performing the run
##################################################
function run!()
    #####
    Ltot_init = get_Ltot()
    Etot_init = get_Etot()
    #####
    print(SEED," |")
    @time integrate_NSTEPS!() # Performing the integration
    flush(stdout) # Forcing the printing
    #####
    Ltot_last = get_Ltot()
    Etot_last = get_Etot()
    err_Ltot = abs((Ltot_last-Ltot_init)/(Ltot_init))
    err_Etot = abs((Etot_last-Etot_init)/(Etot_init))
    println("err_Ltot | ",err_Ltot)
    println("err_Etot | ",err_Etot)
    #####
end
##################################################
run!() # Running the code
