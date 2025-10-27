# The typical value of the frequency is  
# export OPENBLAS_NUM_THREADS=1 <-- ATTENTION, very important to be set, since computing the determinant is the most expensive part
# export OMP_NUM_THREADS=1
# export JULIA_NUM_THREADS=8
# export JULIA_CPU_THREADS=1
#= 
julia --threads 8 run/Run.jl --eps 0.01 --q 0.1 \
--sigma_LOR3 1.0 --J_LOC 1.0 --sigma_LOC 0.2 --T 0.0 --k 1 --nb_int 500
=#
############################################################
include("../src/Main.jl")
############################################################
# It is a good idea to make sure that the two `1D` arrays
# do not have the same size. This can help reducing confusion afterwards
const OMEGA_LEFT = Omega(J_LOC - SIGMA_LOC) # Omega on the left of the action domain. ATTENTION, this number is negative
const OMEGA_CENTER = Omega(J_LOC) # Omega in the centre of the action domain. ATTENTION, this number is negative
const OMEGA_RIGHT = Omega(J_LOC + SIGMA_LOC) # Omega on the right of the action domain. ATTENTION, this number is negative
const OMEGA_IMAG = 10^(-3) * abs(OMEGA_CENTER) # We will make the computation for a tiny POSITIVE imaginary part. ATTENTION, not to forget the absolute value
const RANGE_OMEGA = (OMEGA_RIGHT - OMEGA_LEFT) # Meaningful range in omega. ATTENTION, this number is positive
const TAB_OMEGA_REAL = collect(OMEGA_LEFT - 2 * RANGE_OMEGA : RANGE_OMEGA / 50 : OMEGA_RIGHT + 2 * RANGE_OMEGA) # Table of real(omega). ATTENTION, this range is only good for k=1
const NB_OMEGA = length(TAB_OMEGA_REAL) # Number of real(omega)
const TAB_OMEGA = [ om + OMEGA_IMAG*im for om in TAB_OMEGA_REAL] # Table of omega. ATTENTION, not to forget to add the imaginary part
const TAB_DET = zeros(Complex{Float64},NB_OMEGA) # Table of det
############################################################
println("NB_OMEGA | ",NB_OMEGA)
############################################################
# Function to compute all the needed determinants
############################################################
function TAB_DET!(k::Int64)
    #####
    Threads.@threads for i = 1:NB_OMEGA # Loop over the frequencies
        #####
        omega = TAB_OMEGA[i] # Current omega
        id = Threads.threadid() # Id of the current thread
        tab_E = TAB_E_CALC[id] # Container used for this thread
        val_det = get_detE(k,omega,tab_E) # Computing the determinant
        #####
        TAB_DET[i] = val_det # Filling in the array
        #####
    end
    #####
end
############################################################
# Performing the calculation
############################################################
function run!()
    print("Computing Det | ")
    @time TAB_DET!(K) # Computing the flux for the considered harmonic
    ######
    println("DET [omega=",TAB_OMEGA[100],"] = ",TAB_DET[100])
end
############################################################
run!() # Running