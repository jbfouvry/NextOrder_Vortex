############################################################
# Orbital frequency for the LOR3S profile
# WITHOUT softening
############################################################
function Omega_LOR3S(J::Float64)
    if J <= J_LOR3S # The frequency is zero inside
        return 0.0
    else
        return -0.25*(G*(J-J_LOR3S)*GAMMA_TOT*(J-J_LOR3S+2*SIGMA_LOR3S))/
               (J*pi*(J-J_LOR3S+SIGMA_LOR3S)^(2))
    end
end
############################################################
# Gradient of the orbital frequency, dOmega/dJ,
# for the LOR3S profile -- WITHOUT SOFTENING
############################################################
function OmegaP_LOR3S(J::Float64)
    if J <= J_LOR3S # The gradient of the frequency is zero inside
        return 0.0
    else
        return (G*GAMMA_TOT*((J-J_LOR3S)^(3)+3*(J-J_LOR3S)^(2)*
               SIGMA_LOR3S-2*J_LOR3S*SIGMA_LOR3S^(2)))/(4*J^(2)*
               pi*(J-J_LOR3S+SIGMA_LOR3S)^(3))
    end
end
############################################################
# Function that, given some action J,
# returns, if it exists, an action Jstar
# such that
# (i) Jstar != J
# (ii) Omega(Jstar) = Omega(J)
############################################################
# We compute the action J_FLEX_LOR3S for which
# OmegaP (J_FLEX) = 0
const J_LEFT_LOR3S = J_LOC - SIGMA_LOC
const J_RGHT_LOR3S = J_LOC + SIGMA_LOC
const J_FLEX_LOR3S = bisection(J -> OmegaP_LOR3S(J), J_LEFT_LOR3S, J_RGHT_LOR3S)
############################################################
# We also determine the value of the Omega in [J_LOC-SIGMA_LOC,J_FLEX,J_LOC+SIGMA_LOC]
# We have the ordering
# OMEGA_LEFT, OMEGA_RIGHT > OMEGA_FLEX
const OMEGA_LEFT_LOR3S = Omega_LOR3S(J_LEFT_LOR3S)
const OMEGA_FLEX_LOR3S = Omega_LOR3S(J_FLEX_LOR3S)
const OMEGA_RGHT_LOR3S = Omega_LOR3S(J_RGHT_LOR3S)
############################################################
# In practice, the LOR3S frequency changes of monotonicity
# only once within the interval [J_LOC-SIGMA_LOC, J_LOC+SIGMA_LOC]
# for the action J_FLEX
# So, to find, Jstar, we need
############################################################
function get_Jstar_LOR3S(J::Float64) :: Float64
    #####
    if     (J_LEFT_LOR3S <= J <= J_FLEX_LOR3S) # J is within the left interval
        #####
        # Now, we determine if Omega(J) falls
        # within the frequency domain of the other side
        omega = Omega_LOR3S(J)
        #####
        if (OMEGA_FLEX_LOR3S <= omega <= OMEGA_RGHT_LOR3S) # ATTENTION, to the order
            #####
            Jstar = bisection(Jp -> Omega_LOR3S(Jp)-omega, J_FLEX_LOR3S, J_RGHT_LOR3S)
            #####
            return Jstar # We found a non-trivial resonance
        else
            return NaN # There is no non-trivial resonance
        end
        #####
    elseif (J_FLEX_LOR3S <= J <= J_RGHT_LOR3S) # J is within the right interval
        #####
        # Now, we determine if Omega(J) falls
        # within the frequency domain of the other side
        omega = Omega_LOR3S(J)
        #####
        if (OMEGA_FLEX_LOR3S <= omega <= OMEGA_LEFT_LOR3S) # ATTENTION, to the order
            #####
            Jstar = bisection(Jp -> Omega_LOR3S(Jp)-omega, J_LEFT_LOR3S, J_FLEX_LOR3S)
            #####
            return Jstar # We found a non-trivial resonance
        else
            return NaN # There is no non-trivial resonance
        end
        #####
    end
    #####
    # Otherwise, the code should fail, as we are outside the appropriate domain
    error("get_Jstar_LOR3S: J not within [J_LEFT,J_RGHT]")
    #####
    return NaN
end