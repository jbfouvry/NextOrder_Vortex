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
    # Sampling the massive particles
    # ATTENTION, to the mass of the individual particles
    for i=1:NMASS
        gamma = Q * GAMMATOT / NMASS # Initial circulation of the particle. ATTENTION, to the active fraction Q
        x,y = init_1part() # Drawing one particle
        # Filling in the arrays
        TAB_GAMMA[i] = gamma
        TAB_X[i]     = x
        TAB_Y[i]     = y
    end
    #####
    # Samping the test particles
    for i=(NMASS+1):NPART
        gamma = 0.0 # Mass of the test particles
        phi = 2 * pi * rand() # Random phase of the particle
        r = sqrt(2 * JTEST) # Initial radius of the test particle
        s, c = sincos(phi)
        x, y = r * c, r * s # Initial position of the test particle
        #####
        # Filling in the arrays
        TAB_GAMMA[i] = gamma
        TAB_X[i] = x
        TAB_Y[i] = y
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