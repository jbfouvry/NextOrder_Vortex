##################################################
# Function to recenter the simulation so that
# rtot = SUM_i {x_i,y_i} = {0,0} initially
##################################################
function init_center!()
    rtot_x, rtot_y = get_rtot() # Computing the initial invariant r_tot = SUM_i[\gamma_i {x_i,y_i}]
    # Rescaling the invariant so that it has the dimension of a length
    # ATTENTION, not to forget the active fraction
    # @IMPROVE -- This would not work with only test particles
    rtot_x /= (Q*GAMMATOT)
    rtot_y /= (Q*GAMMATOT)
    #####
    for i=1:NPART
        TAB_X[i] -= rtot_x
        TAB_Y[i] -= rtot_y
    end
end
##################################################
# Picking the sampling function
##################################################
const init_1part = init_1part_LOC
##################################################
# Sampling the particles
##################################################
function init!()
    for i=1:NPART
        gamma = Q*GAMMATOT/NPART # Initial circulation of the particle. ATTENTION, to the active fraction Q
        x,y = init_1part() # Drawing one particle
        # Filling in the arrays
        TAB_GAMMA[i] = gamma
        TAB_X[i]     = x
        TAB_Y[i]     = y
    end
    #####
    # We recentre the initial sampling
    if (RECENTER == "true")
        init_center!() # Recentering the initial conditions
    end
end
##################################################
init!() # Initialising the particles
##################################################
# To accelerate the code, we use the fact that all the particles
# have the same mass, in TAB_DOT!()
##################################################
const GAMMA = TAB_GAMMA[1] # Mass of the first particle
##################################################
function check_GAMMA!()
    #####
    # Checking that the mass of the particle 1 is non-zero
    if (GAMMA == 0.0)
        error("check_GAMMA! -- The mass of the test particle is 0")
    end
    #####
    # Checking that all the particles have the same mass
    for n=1:NPART
        if (TAB_GAMMA[n] != GAMMA)
            error("check_GAMMA! -- The particles do not have al the same mass")
        end
    end
end
##################################################
check_GAMMA!() # Checking that all the particles have the mass
