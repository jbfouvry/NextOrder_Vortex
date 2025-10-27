############################################################
# Wrapped function that computes the rescaled flux
# for a given action, while computing the coupling coefficient
# up to a certain kmax
############################################################
function FluxRESCALED(kmax::Int64,J::Float64)  :: Float64
    #####
    Jstar = get_Jstar(J) # Solving for the resonance condition
    #####
    # Making sure that this is not a NaN,
    # i.e. there is a non-trivial resonance
    if (Jstar == Jstar) # Making sure that this is not a NaN,
        #####
        return 2 * pi^(2) * UtotSQ(kmax,J,Jstar) *
               (DF(Jstar) * DFP(J) - DF(J) * DFP(Jstar)) /
               (abs(OmegaP(Jstar)))  # ATTENTION, not to forget the absolute value
        #####
        return res
    #####
    # Otherwise, there is no non-trivial resonance
    # and the flux is 0.0
    else
        return 0.0
    end
end
