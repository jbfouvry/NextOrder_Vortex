############################################################
# Mean DF for the LOC profile
# ATTENTION, not to forget the active fraction
# In order to improve the code's accuracy,
# we also introduce the function (dF/dJ)/F
# This will help stabilising the computation of the
# crossed term (k.d/dJ)F_3 in the kinetic equation
# ATTENTION, not to forget the active fraction
############################################################
function DF_LOC(J::TYPE)
    #####
    x = (J - J_LOC) / SIGMA_LOC # Rescaled distance to the center
    #####
    if (x <= -1) || (x >= 1) # We are outside the definition domain
        return TYPE(0)
    end
    #####
    pref = Q * 15 * GAMMA_TOT / (32 * PI * SIGMA_LOC) # Prefactor. ATTENTION, to the active fraction
    #####
    return pref * (x - 1)^(2) * (x + 1)^(2)
end
#####
function DFP_LOC(J::TYPE)
    #####
    x = (J - J_LOC) / SIGMA_LOC # Rescaled distance to the center
    #####
    if (x <= -1) || (x >= 1) # We are outside the definition domain
        return TYPE(0)
    end
    #####
    pref = Q * 15 * GAMMA_TOT / (8 * PI * SIGMA_LOC^(2)) # Prefactor. ATTENTION, to the active fraction
    #####
    return pref * x * (x - 1) * (x + 1)
end
#####
# Ratio F'/F. ATTENTION, this diverges at the edge of the definition domain
#####
function DFPoverDF_LOC(J::TYPE)
    #####
    x = (J - J_LOC) / SIGMA_LOC # Rescaled distance to the center
    #####
    if (x <= -1) || (x >= 1) # We are outside the definition domain
        return TYPE(0)
    end
    #####
    pref = 4 / SIGMA_LOC # Prefactor. ATTENTION, this does NOT contain the active fraction
    #####
    return pref * x / (x + 1) / (x - 1)
end