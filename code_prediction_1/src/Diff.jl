############################################################
# Wrapped function that computes the diffusion coefficient
# for a given action, while computing the coupling coefficient
# ATTENTION, this is the rescaled diffusion coefficient,
# i.e. D / gamma
############################################################
function calc_Diff(J::Float64) :: Float64
    #####
    # Adding the local contribution Jstar = J
    res = 4 * pi^(2) * UtotSQ(KMAX,J,J) * DF(J) / (abs(OmegaP(J))) # ATTENTION, not to forget the absolute value
    #####
    # Now, we determine if there exists a non-local resonance condition
    Jstar = get_Jstar(J) # Solving for the resonance condition
    #####
    # Making sure that this is not a NaN,
    # i.e. there is a non-trivial resonance
    if (Jstar == Jstar) # Making sure that this is not a NaN,
        #####
        res += 4 * pi^(2) * UtotSQ(KMAX,J,Jstar) * DF(Jstar) / (abs(OmegaP(Jstar)))  # ATTENTION, not to forget the absolute value
        #####
    end
    #####
    return res
end
############################################################
# Wrapped function that computes the Landau diffusion coefficient
# on the grid TAB_J_INT
############################################################
const TAB_DIFF = zeros(Float64,NB_INT) # Container for the Landau D
#####
function TAB_DIFF!()
    for i_J=1:NB_INT
        J_loc = TAB_J_INT[i_J]
        D_loc = calc_Diff(J_loc)
        TAB_DIFF[i_J] = D_loc
    end
end


