############################################################
# Mean potential and orbital frequency for the LOR3 profile
# WITHOUT softening
############################################################
const H_0_LOR3 = (G*GAMMATOT)/(4.0*pi) # Scale for the Hamiltonian
############################################################
function U_LOR3(J::Float64)
    x = J / SIGMA_LOR3 # Rescaled action
    return H_0_LOR3*(1.0-(1.0+x)*log(2.0*(J+SIGMA_LOR3)))/(1.0+x)
end
############################################################
const OMEGA_0_LOR3 = - (G*GAMMATOT)/(2.0*pi*SIGMA_LOR3) # Scale for the frequency. ATTENTION, to the minus sign
############################################################
function _g_LOR3(x::Float64)
    return (1.0+0.5*x)/((1.0+x)^(2))
end
############################################################
function Omega_LOR3(J::Float64)
    return OMEGA_0_LOR3*_g_LOR3(J/SIGMA_LOR3)
end