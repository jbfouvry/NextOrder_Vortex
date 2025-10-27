############################################################
# Wrapped function that computes the diffusion coefficient
# by solving the self-consistent broadened equation it satisfies
############################################################
const DELTA_J_INT = (2 * SIGMA_LOC) / NB_INT # Step distance to perform integrals
const TAB_J_INT = [J_LOC - SIGMA_LOC + (i_J - 0.5) * DELTA_J_INT for i_J = 1:NB_INT] # Table of the nodes in u
const TAB_DIFF_SELF = zeros(Float64,NB_INT) # Container for the self-consistent D
const TAB_DIFF_SELF_NEW = zeros(Float64,NB_INT) # Temporary container for D
############################################################
# Function that initialises the self-consistent diffusion coefficient
# @TODO -- Should be initialised with the correct scaling, at least
############################################################
function TAB_DIFF_SELF_INIT!()
    for i_J = 1:NB_INT
        D_init = Q * GAMMA_PART # Initial (lazy) guess for the diffusion coefficient.
        TAB_DIFF_SELF[i_J] = D_init # Filling in the array
    end
end
############################################################
# Returns the diffusion coefficient D_thetatheta,
# from the local value of D_JJ, at the location J,
# assuming an isotropic diffusion
############################################################
function get_Dthetatheta(Diff,J)
    return Diff / (4 * J^2)
end
############################################################
# Returns the local regularising time
# for a given harmonic, k, and diffusion coefficient, D
############################################################
# const FACTOR_FIDUCIAL_TREG = 0.0001
const FACTOR_FIDUCIAL_TREG = 1.0
#####
function get_Treg_SELF(k,Diff,J)
    D_thetatheta = get_Dthetatheta(Diff,J) # Local value of the diffusion coefficient in theta
    return FACTOR_FIDUCIAL_TREG / (D_thetatheta * k^(2))
end
############################################################
# Returns the value of the broadened Dirac delta
# for a given regularising time, T, and frequency, omega
############################################################
function get_delta_SELF(T,omega)
    return (1 / pi) * T / (1 + (omega * T)^(2))
end
############################################################
# Computes the new diffusion coefficient at a location J,
# given the 'background' diffusion coefficient stored in tab_D
# @IMPROVE -- Maybe better to perform the sum over k inside
############################################################
function calc_Diff_self(J,tab_Diff)
    #####
    res = 0.0 # Initialising the result
    #####
    omega = Omega(J) # Value of the frequency in J
    #####
    for k=1:KMAX # Loop over the harmonic indices
        #####
        res_loc = 0.0
        #####
        for i_J = 1:NB_INT # Performing the integral
            #####
            J1 = TAB_J_INT[i_J] # Local u
            F1 = DF(J1) # Local value of the DF
            omega1 = Omega(J1) # Local value of Omega
            k_DeltaOmega = k * (omega - omega1) # Value of the frequency difference
            U_sq = UPotk(k,J,J1)^(2) # Coupling coefficient squared
            #####
            Diff1 = tab_Diff[i_J] # Local value of the 'background' D
            T1 = get_Treg_SELF(k,Diff1,J1) # Local regularising time
            #####
            delta_broad = get_delta_SELF(T1,k_DeltaOmega) # Evaluating the broadened Dirac
            #####
            res_loc += delta_broad * F1 * U_sq # Adding the contribution
            ##### 
        end
        #####
        # Multiplying by the prefactor of the sum over k
        pref_loc = k^(2) * DELTA_J_INT # Prefactor in the sum over k
        res_loc *= pref_loc
        #####
        res += res_loc
    end
    #####
    # Global prefactor of the expression
    # ATTENTION, to the prefactor "8"
    pref = 8 * pi^(2) * GAMMA_PART
    res *= pref
    #####
    return res # Output
end
############################################################
# Performing one iteration of the fix-point search
# by updating TAB_DIFF_SELF_NEW
# ATTENTION, we have added parallelisation
############################################################
function TAB_DIFF_SELF_NEW!()
    Threads.@threads for i_J = 1:NB_INT
        #####
        J = TAB_J_INT[i_J] # Location in u for which we want to compute D
        #####
        # Predicting a new coefficient Diff, at the location u,
        # with the 'background' Diff stored in TAB_DIFF_SELF
        Diff_new = calc_Diff_self(J,TAB_DIFF_SELF)
        #####
        TAB_DIFF_SELF_NEW[i_J] = Diff_new # Updating the container
    end
end
############################################################
# Returns the (relative) distance between TAB_D_SELF and TAB_D_SELF_NEW
# It is defined as the maximum relative distance, coordinate-wise,
# using TAB_D_SELF as the truth
############################################################
function get_dist_SELF()
    #####
    dist = - Inf # Initialising the distance
    #####
    for i_J = 1:NB_INT # Running over all the coordinates
        Diff = TAB_DIFF_SELF[i_J] # Old D
        Diff_new = TAB_DIFF_SELF_NEW[i_J] # New D
        #####
        err = abs((Diff_new - Diff) / Diff) # Relative distance, using the old D as the truth
        #####
        if (err > dist) # We have found a worse distance
            dist = err
        end
    end
    #####
    return dist # Output
end
############################################################
# Utility function that copies the content of TAB_DIFF_SELF_NEW
# into TAB_DIFF_SELF
############################################################
function TAB_DIFF_SELF_update!()
    for i_J = 1:NB_INT
        TAB_DIFF_SELF[i_J] = TAB_DIFF_SELF_NEW[i_J] # Updating the content of the array
    end
end
############################################################
# Wrapped function that performs iteration
# for the SELF diffusion coefficient, up to convergence
# Here, the stopping criteria is an relative distance between iterations
############################################################
function TAB_DIFF_SELF!()
    #####
    TAB_DIFF_SELF_INIT!() # Initialisation
    #####
    iter = 0 # Number of iterations
    dist = + Inf # Distance between iterations
    #####
    # Performing iterations, up to convergence
    while (dist > DIST_SELF) && (iter < MAX_ITER_SELF)
        #####
        TAB_DIFF_SELF_NEW!() # Computing a new diffusion coefficient
        iter += 1 # We have performed one more iteration
        #####
        dist = get_dist_SELF() # Getting the distance between the two diffusion coefficients
        println("iter | ",iter," | ",dist)
        #####
        TAB_DIFF_SELF_update!() # Copying the new array into the old one
        #####
    end
    #####
    println("DIFF_SELF | iter | ",iter," | ",dist)
    #####
end
