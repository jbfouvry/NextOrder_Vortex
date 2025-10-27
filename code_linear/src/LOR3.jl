############################################################
# Orbital frequency and gradients
############################################################
const OMEGA_0_LOR3 = - (G * GAMMA_TOT) / (2 * pi * SIGMA_LOR3) # Scale for the frequency
#####
function _g_LOR3(x::Float64)
    return (1 + (x / 2)) / ((1 + x)^(2))
end
#####
function Omega_LOR3(J::Float64)
    return OMEGA_0_LOR3 * _g_LOR3(J / SIGMA_LOR3)
end