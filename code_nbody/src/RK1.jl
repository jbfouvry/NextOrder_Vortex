##################################################
# Integrates for one timestep using RK1
##################################################
function integrate_DT_RK1!()
    #####
    TAB_DOT!() # Computing the velocities
    #####
    drift!(TAB_X,TAB_Y,
           TAB_X,TAB_Y,
           TAB_XDOT,TAB_YDOT,
           DT) # Full drift
    #####
    TIME[1] += 1 # Updating the time
end