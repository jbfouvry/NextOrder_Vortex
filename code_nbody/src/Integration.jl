const TIME = zeros(Int64, 1) # Global time of the simulation
TIME[1] = 0 # Initialising the time
##################################################
# Copying the initial location of the particles
# of the content of (tab_x_0,tab_y_0) into
# (tab_x_1,tab_y_1)
##################################################
function copy_tab!(tab_x_0,tab_y_0,tab_x_1,tab_y_1)
    for i=1:NPART
        tab_x_1[i] = tab_x_0[i]
        tab_y_1[i] = tab_y_0[i]
    end
end
##################################################
# Combine the content of arrays
# with given prefactors
##################################################
function combine_tab!(list_coeff,
                      list_tab_X,
                      list_tab_Y,
                      tab_out_X,
                      tab_out_Y)
    #####
    nb_terms = length(list_coeff) # Number of tab that have to be combined
    ######
    # Initialising the arrays with 0.0
    fill!(tab_out_X,0.0)
    fill!(tab_out_Y,0.0)
    ######
    for i=1:nb_terms # Loop over the terms
        coeff = list_coeff[i] # Current coeff
        #####
        tab_loc_X = list_tab_X[i] # Current X array
        tab_loc_Y = list_tab_Y[i] # Current Y array
        for n=1:NPART
            tab_out_X[n] += coeff*tab_loc_X[n]
            tab_out_Y[n] += coeff*tab_loc_Y[n]
        end
    end
end
##################################################
# Drifting all the particles
# starting from (tab_x_0,tab_y_0)
# with the velocities (tab_xdot,tab_ydot)
# filling in (tab_x_1,tab_y_1)
# for a duration dt
##################################################
function drift!(tab_x_0,tab_y_0,
                tab_x_1,tab_y_1,
                tab_xdot,tab_ydot,
                dt::Float64)
    #####
    for i=1:NPART
        tab_x_1[i] = tab_x_0[i] + dt*tab_xdot[i]
        tab_y_1[i] = tab_y_0[i] + dt*tab_ydot[i]
    end
end
##################################################
# Performs an integration for NSTEPS steps
##################################################
function integrate_NSTEPS!()
    for istep=1:NSTEPS
        integrate_DT!()
    end
end
##################################################
# Returns the current time
##################################################
function get_time()
    return TIME[1]*DT
end
##################################################
if     (SCHEME == "RK1")
    include("RK1.jl")
    const integrate_DT! = integrate_DT_RK1!
elseif (SCHEME == "RK2")
    include("RK2.jl")
    const integrate_DT! = integrate_DT_RK2!
elseif (SCHEME == "RK4")
    include("RK4.jl")
    const integrate_DT! = integrate_DT_RK4!
elseif (SCHEME == "RK6")
    include("RK6.jl")
    const integrate_DT! = integrate_DT_RK6!
elseif (SCHEME == "RK9")
    include("RK9.jl")
    const integrate_DT! = integrate_DT_RK9!
elseif (SCHEME == "RK14")
    include("RK14.jl")
    const integrate_DT! = integrate_DT_RK14!
end