# Temporary containers for the integrations
const TAB_X_1 = zeros(Float64,NPART)
const TAB_Y_1 = zeros(Float64,NPART)
#####
const TAB_XDOT_1 = zeros(Float64,NPART)
const TAB_YDOT_1 = zeros(Float64,NPART)
const TAB_XDOT_2 = zeros(Float64,NPART)
const TAB_YDOT_2 = zeros(Float64,NPART)
const TAB_XDOT_3 = zeros(Float64,NPART)
const TAB_YDOT_3 = zeros(Float64,NPART)
const TAB_XDOT_4 = zeros(Float64,NPART)
const TAB_YDOT_4 = zeros(Float64,NPART)
##################################################
# Integrates for one timestep using RK4
##################################################
function integrate_DT_RK4!()
    #####
    copy_tab!(TAB_X,TAB_Y,TAB_X_1,TAB_Y_1) # Copying the initial location
    #####
    TAB_DOT!()
    copy_tab!(TAB_XDOT,TAB_YDOT,TAB_XDOT_1,TAB_YDOT_1)
    #####
    drift!(TAB_X_1,TAB_Y_1,
           TAB_X,TAB_Y,
           TAB_XDOT,TAB_YDOT,
           0.5*DT)
    #####
    TAB_DOT!()
    copy_tab!(TAB_XDOT,TAB_YDOT,TAB_XDOT_2,TAB_YDOT_2)
    #####
    drift!(TAB_X_1,TAB_Y_1,
           TAB_X,TAB_Y,
           TAB_XDOT,TAB_YDOT,
           0.5*DT)
    #####
    TAB_DOT!()
    copy_tab!(TAB_XDOT,TAB_YDOT,TAB_XDOT_3,TAB_YDOT_3)
    #####
    drift!(TAB_X_1,TAB_Y_1,
           TAB_X,TAB_Y,
           TAB_XDOT,TAB_YDOT,
           DT)
    #####
    TAB_DOT!()
    copy_tab!(TAB_XDOT,TAB_YDOT,TAB_XDOT_4,TAB_YDOT_4)
    #####
    combine_tab!((1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0),
              (TAB_XDOT_1, TAB_XDOT_2, TAB_XDOT_3, TAB_XDOT_4),
              (TAB_YDOT_1, TAB_YDOT_2, TAB_YDOT_3, TAB_YDOT_4),
               TAB_XDOT,
               TAB_YDOT)
    #####
    drift!(TAB_X_1,TAB_Y_1,
           TAB_X,TAB_Y,
           TAB_XDOT,TAB_YDOT,
           DT)
    #####
    TIME[1] += 1 # Updating the time
end