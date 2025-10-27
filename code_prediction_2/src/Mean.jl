############################################################
# Functions specific to the considered initial conditions
############################################################
# Definition of the background
############################################################
include("LOR3.jl")
#####
const OMEGA_0 = OMEGA_0_LOR3 # Frequency scale
const _g = _g_LOR3 # Frequency function
const _dgdx = _dgdx_LOR3 # Gradient of the frequency
const _diffg = _diffg_LOR3 # Difference of frequency
const Jres = Jres_LOR3 # Resonant action
############################################################
# Definition of the active DF
############################################################
include("LOC.jl")
#####
const DF = DF_LOC # Active DF
const DFP = DFP_LOC # Gradient of the DF
const DFPoverDF = DFPoverDF_LOC # Ratio F'/F
############################################################
# Orbital frequency and gradients
# ATTENTION, this is rescaled wrt SIGMA_LOR3
############################################################
function Omega(J::TYPE)
    return OMEGA_0 * _g(J / SIGMA_LOR3)
end
#####
function OmegaP(J::TYPE)
    return (OMEGA_0 / SIGMA_LOR3) * _dgdx(J / SIGMA_LOR3)
end
#####
# Computes Omega[J]-Omega[J1]
# hopefully in a numerically more stable fashion
function diffOmega(J::TYPE,J1::TYPE)
    return OMEGA_0 * _diffg(J / SIGMA_LOR3, J1 / SIGMA_LOR3)
end
############################################################
const TDYN = abs(2 * pi / Omega(J_LOC)) # Dynamical time