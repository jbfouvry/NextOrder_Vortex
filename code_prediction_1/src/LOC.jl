############################################################
# Mean DF for the LOC profile
# ATTENTION, not to forget the active fraction
# In order to improve the code's accuracy,
# we also introduce the function (dF/dJ)/F
# This will help stabilising the computation of the
# crossed term (k.d/dJ)F_3 in the kinetic equation
# ATTENTION, not to forget the active fraction
############################################################
function DF_LOC(J::Float64)
    #####
    x = (J - J_LOC) / SIGMA_LOC # Rescaled distance to the center
    #####
    if (x <= -1) || (x >= 1) # We are outside the definition domain
        return 0.0
    end
    #####
    pref = Q * 15 * GAMMA_TOT / (32 * pi * SIGMA_LOC) # Prefactor. ATTENTION, to the active fraction
    #####
    return pref * (x - 1)^(2) * (x + 1)^(2)
end
#####
function DFP_LOC(J::Float64)
    #####
    x = (J - J_LOC) / SIGMA_LOC # Rescaled distance to the center
    #####
    if (x <= -1) || (x >= 1) # We are outside the definition domain
        return 0.0
    end
    #####
    pref = Q * 15 * GAMMA_TOT / (8 * pi * SIGMA_LOC^(2)) # Prefactor. ATTENTION, to the active fraction
    #####
    return pref * x * (x - 1) * (x + 1)
end