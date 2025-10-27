############################################################
# Gradient of the Mean DF for the LOC profile
# ATTENTION, not to forget the active fraction
############################################################
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