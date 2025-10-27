##################################################
# Function to compute the velocities of all the particles
# ATTENTION -- Here, we assume that the particles can have different mass
# @IMPROVE -- Here, we sum the force from the test particles
#             even if their mass is 0.0
##################################################
@fastmath function TAB_DOT!()
    #####
    for i=1:NPART
        @inbounds x_i, y_i = TAB_X[i], TAB_Y[i]
        J_i = 0.5*(x_i^(2) + y_i^(2)) # Action of the particle
        Om_i = Omega(J_i) # Getting the externally-imposed frequency of the particle
        #####
        # ATTENTION, we do not use "+=", but make the first fill of the arrays
        @inbounds TAB_XDOT[i] =  Om_i*y_i
        @inbounds TAB_YDOT[i] = -Om_i*x_i
    end
    #####
    for i=1:NPART
        @inbounds gamma_i, x_i, y_i = TAB_GAMMA[i], TAB_X[i], TAB_Y[i] # Circulation and position of particle i
        acc_x, acc_y = 0.0, 0.0 # Accumulators for the force on particle i
        for j=(i+1):NPART
            @inbounds gamma_j, x_j, y_j = TAB_GAMMA[j], TAB_X[j], TAB_Y[j] # Circulation and position of particle j
            #####
            inv_norm = (G/(2 * pi))/((x_i - x_j)^(2) + (y_i - y_j)^(2) + EPS^(2)) # Inverse of the distance squared between the particles. ATTENTION, this is multiplied by the appropriate prefactor. ATTENTION, we have added the softening
            #####
            dxdt = -(y_i - y_j)*inv_norm # Up to gamma_j, this is dxi/dt. ATTENTION, to the minus sign
            dydt =  (x_i - x_j)*inv_norm # Up to gamma_j, this is dyi/dt
            #####
            acc_x += gamma_j*dxdt # Contribution to the force on particle i
            acc_y += gamma_j*dydt # Contribution to the force on particle i
            #####
            @inbounds TAB_XDOT[j] -= gamma_i*dxdt # ATTENTION, to the minus sign
            @inbounds TAB_YDOT[j] -= gamma_i*dydt # ATTENTION, to the minus sign
        end
        #####
        # Adding to the force on particle i
        @inbounds TAB_XDOT[i] += acc_x
        @inbounds TAB_YDOT[i] += acc_y
    end
end
