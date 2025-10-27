############################################################
# Wrapped function that computes the RESCALED diffusion coefficient
# for a given action
# ATTENTION, we regularise the DIRAC *after* performing the sum over k
# ATTENTION, this is the rescaled diffusion coefficient,
# i.e. D / gamma
############################################################
function calc_Diff_REG(kmax::Int64,T::Float64,J::Float64)
    #####
    res = 0.0 # Initialising the result
    #####
    omega = Omega(J) # Frequency
    #####
    for i_J = 1:NB_INT
        #####
        Jp = TAB_J_INT[i_J] # Current action
        #####
        omegap = Omega(Jp) # Current frequency
        #####
        res += UtotSQ(kmax,J,Jp) * DF(Jp) *
               kernel_REG_LRZ(T, omega - omegap) # ATTENTION, not yet multiplied by the prefactor
        #####
    end
    #####
    res *= DELTA_J_INT * 4 * pi^(2) # ATTENTION, need to multiply by DeltaJ. ATTENTION, to the factor 2
    #####
    return res
end
############################################################
# Wrapped function that computes the regularised diffusion coefficient
# on the grid TAB_J_INT
############################################################
const TAB_DIFF_REG = zeros(Float64,NB_INT) # Container for the Landau D
#####
function TAB_DIFF_REG!()
    for i_J=1:NB_INT
        J_loc = TAB_J_INT[i_J]
        D_loc = calc_Diff_REG(KMAX,TREG,J_loc)
        TAB_DIFF_REG[i_J] = D_loc
    end
end
