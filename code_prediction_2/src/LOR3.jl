############################################################
# Orbital frequency and gradients
############################################################
const OMEGA_0_LOR3 = - (G * GAMMA_TOT) / (2 * PI * SIGMA_LOR3) # Scale for the frequency
#####
function _g_LOR3(x::TYPE)
    return (1 + (x / 2)) / ((1 + x)^(2))
end
#####
function _dgdx_LOR3(x::TYPE)
    return - (3 + x) / (2 * (1 + x)^(3))
end
#####
# Computes the difference between the frequencies
# This provides a slightly more stable numerical evaluation
function _diffg_LOR3(x::TYPE,xp::TYPE)
    return - (x-xp) * (3 + 2 * (x+xp) + x*xp) / (2 * (1 + x)^(2)*(1 + xp)^(2))
end
############################################################
# Function that returns, if it exists
# the resonant action Jres such that Omega[Jres]=omega
# If the resonance cannot be solved, we return NaN
############################################################
function Jres_LOR3(omega::TYPE)
    #####
    om_res = omega / OMEGA_0_LOR3 # Rescaled resonant frequency
    # If it is in ]0,1], then we can solve the resonance condition
    # ATTENTION, we do not allow for resonance exactly at zero frequency
    # Those can sometimes occur because of the simple rational fraction expression
    # of the LOR3 frequency
    if (0 < om_res <= 1) # We can hope to solve for the resonance condition
        x_res = (1 - 4 * om_res + sqrt(1 + 8 * om_res)) / (4 * om_res) # Resonant dimensionless action, x=J/SIGMA_LOR3
        J_res = SIGMA_LOR3 * x_res # Resonant action
        return J_res # Output
    else # There is no resonant action
        return NaN # We return a NaN
    end
end