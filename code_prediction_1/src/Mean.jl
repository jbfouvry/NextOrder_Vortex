############################################################
# Definition of the active DF
############################################################
include("LOC.jl")
#####
const DF = DF_LOC # Active DF
const DFP = DFP_LOC # Gradient of the DF
############################################################
# Orbital frequency and gradients,
# and resonance condition solver
############################################################
include("LOR3S.jl")
const Omega = Omega_LOR3S
const OmegaP = OmegaP_LOR3S
const get_Jstar = get_Jstar_LOR3S
############################################################
const TDYN = abs(2 * pi / Omega(J_LOC)) # Dynamical time
const TREG = TREG_OVER_TDYN * TDYN # Regularising time in physical units