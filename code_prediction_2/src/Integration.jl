############################################################
# Wrapped function that computes the rescaled flux
# for a given fundamental resonance {k,kp}
# and a given action J
# The integral has the form
# INT[IntFluxRescaled[k,kp,J,J1],{J1}]
# This integral is difficult to compute because it involves
# a resonant denominator that must be computed following
# Landau's prescription
# We perform the integration in log-space to always have
# an appropriate resolution
# We always sum in a symmetric fashion, i.e. following the
# prescription from Cauchy's principal value
# In practice, we write
# I = \int_{Jm}^{Jp} dJ1 F[J,J1]
# with Jm = J_LOC - SIGMA_LOC
#      Jp = J_LOC + SIGMA_LOC
# We rewrite this integral as
# I = \int_{0}^{Max[Dm,Dp]} dD (f[J,J-D] * Theta[D <= Dm]+
#                               f[J,J+D] * Theta[D <= Dp])
# with Dm = J  - Jm
#      Dp = Jp - J
# In practice, we evaluate the integral using a midpoint rule
# with a fixed step given by DI = Max[Dm,Dp] / NB_INT
# ATTENTION, it is essential not to use too many sampling nodes
# otherwise the cancellation in Cauchy's symmetrisation
############################################################
function FluxRESCALED(k::Int64,kp::Int64,J::TYPE)
    #####
    res = TYPE(0) # Initialising the result
    #####
    # Bounds of the action domain
    Jm = J_LOC - SIGMA_LOC
    Jp = J_LOC + SIGMA_LOC
    #####
    # Moving to D-space ["D" stands for "Delta"]
    Dm = J - Jm
    Dp = Jp - J
    Dmax = max(Dm,Dp) # Upper bound for the integration
    #####
    # For safety, we check that these numbers are real
    @assert (Dm > 0) "FluxRESCALED -- must have Dm > 0"
    @assert (Dp > 0) "FluxRESCALED -- must have Dp > 0"
    #####
    # Step distance considered
    DI = Dmax / NB_INT
    #####
    for i=1:NB_INT # Loop over the sampling points
        #####
        D = TYPE(0) + DI * (i - (1//2)) # Current value of D=Delta. ATTENTION, to the '1//2', we use midpoint
        #####
        # Minus value
        J1m = J - D
        if (Jm < J1m < Jp) # We are within the DF domain
            vm = IntFluxRESCALED(k,kp,J,J1m)
        else
            vm = TYPE(0)
        end
        #####
        # Plus value
        J1p = J + D
        if (Jm < J1p < Jp)
            vp = IntFluxRESCALED(k,kp,J,J1p)
        else
            vp = TYPE(0)
        end
        #####
        v = vm + vp # Summing the left/rght contributions
        #####
        res += v # Adding the contributions to the flux
        #####
    end
    #####
    res *= DI # Rescaling by the step size
    #####
    return res # Output
end
