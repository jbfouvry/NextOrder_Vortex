##################################################
const TAB_GAMMA = zeros(Float64,NPART) # Individual masses
const TAB_X     = zeros(Float64,NPART) # X of the particles
const TAB_Y     = zeros(Float64,NPART) # Y of the particles
const TAB_XDOT = zeros(Float64,NPART) # Rates of change, dx/dt
const TAB_YDOT = zeros(Float64,NPART) # Rates of change, dy/dt
############################################################
# Definition of the BACKGROUND potential and frequency
############################################################
if BACKGROUND == "LOR3"
    include("LOR3.jl")
    const U_back = U_LOR3
    const Omega_back = Omega_LOR3
elseif (BACKGROUND == "LOR3S")
    include("LOR3S.jl")
    const U_back = U_LOR3S
    const Omega_back = Omega_LOR3S
end
############################################################
# Mean external potential and external frequency
# Here, we remove the (average) mean potential generated with softening,
# namely LOC_soft,
# and add the (unsoftened) mean potential we want to reach,
# namely LOR3
############################################################
function Uext(J::Float64)
    return U_back(J) - Q * U_LOC_soft(J)
end
############################################################
function Omega(J::Float64)
    return Omega_back(J) - Q * Omega_LOC_soft(J)
end
############################################################
# Computes the total momentum of the system
##################################################
function get_Ltot()
    #####
    Ltot = 0.0 # Initialising the counter
    #####
    for i=1:NPART # Loop over the particles
        x_i, y_i = TAB_X[i], TAB_Y[i] # Position of the particle
        r_i_SQ = x_i^(2) + y_i^(2) # Radius squared of the particle
        Ltot += TAB_GAMMA[i]*r_i_SQ # Adding the angular momentum from particle i
    end
    #####
    return Ltot # Output
end
##################################################
# Computes the total momentum of the system
##################################################
function get_rtot()
    #####
    x_tot, y_tot = 0.0, 0.0 # Initialising the counter
    #####
    for i=1:NPART # Loop over the particles
        gamma_i, x_i, y_i = TAB_GAMMA[i], TAB_X[i], TAB_Y[i] # Circulation and position of the particle
        # Adding the contribution from particle i
        x_tot += gamma_i*x_i
        y_tot += gamma_i*y_i
    end
    #####
    return x_tot, y_tot # Output
end
##################################################
# Computes the total energy of the system
# The specific pairwise interaction potential is
# U = - G/(2pi)ln(|r-r'|)
#   = - G/(4pi)ln(|r-r'|^2)
##################################################
function get_Etot()
    #####
    # External potential
    E_ext = 0.0
    #####
    for i=1:NPART
        gamma_i, x_i, y_i = TAB_GAMMA[i], TAB_X[i], TAB_Y[i]
        J_i = 0.5*(x_i^(2) + y_i^(2)) # Action of the particle
        E_ext += gamma_i * Uext(J_i)
    end
    #####
    # Self-consistent potential
    E_self = 0.0
    #####
    for i=1:NPART
        gamma_i, x_i, y_i = TAB_GAMMA[i], TAB_X[i], TAB_Y[i]
        for j=(i+1):NPART
            gamma_j, x_j, y_j = TAB_GAMMA[j], TAB_X[j], TAB_Y[j]
            #####
            r_ij_SQ = (x_i-x_j)^(2) + (y_i-y_j)^(2) + EPS^(2) # Squared distance between the particles. ATTENTION, this is softened
            E_self -= gamma_i*gamma_j*log(r_ij_SQ) # Adding the contribution. ATTENTION, to the minus sign. ATTENTION, the prefactor is not yet added
        end
    end
    #####
    E_self *= G/(4.0*pi) # Multiplying by the overall prefactor. ATTENTION, here we use 1/(4pi) because we sum distance-squared
    #####
    E_tot = E_ext + E_self # Adding the external and self-gravitating components
    #####
    return E_tot # Output
end
