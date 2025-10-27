const J_MIN = J_LOC - SIGMA_LOC # Minimum J considered for the DF
const J_MAX = J_LOC + SIGMA_LOC # Maximum J considered for the DF
const DELTA_J = (J_MAX - J_MIN) / (NB_INT) # Step distance considered
############################################################
# Background profile
if BACKGROUND == "LOR3"
    include("LOR3.jl")
    const Omega = Omega_LOR3
elseif (BACKGROUND == "LOR3S")
    include("LOR3S.jl")
    const Omega = Omega_LOR3S
end
############################################################
# Active DF
include("LOC.jl")
const DFP = DFP_LOC # Gradient of the DF
############################################################
const TAB_J = [J_MIN + DELTA_J * (i - (1//2)) for i=1:NB_INT] # Table of the discrete J considered 
