# Temporary containers for the integrations
const TAB_X_1 = zeros(Float64,NPART)
const TAB_Y_1 = zeros(Float64,NPART)
##################################################
# Integrates for one timestep using RK2
##################################################
function integrate_DT_RK2!()
    #####
    copy_tab!(TAB_X,TAB_Y,TAB_X_1,TAB_Y_1) # Copying the initial location
    #####
    TAB_DOT!() # Computing the velocities
    #####
    drift!(TAB_X_1,TAB_Y_1,
           TAB_X,TAB_Y,
           TAB_XDOT,TAB_YDOT,
           0.5*DT) # Drifting for half a timestep
    #####
    TAB_DOT!()
    #####
    drift!(TAB_X_1,TAB_Y_1,
           TAB_X,TAB_Y,
           TAB_XDOT,TAB_YDOT,
           DT) # Full drift
    #####
    TIME[1] += 1 # Updating the time
end