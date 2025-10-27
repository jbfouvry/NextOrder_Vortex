const TAB_E_CALC = [zeros(Complex{Float64},NB_INT,NB_INT) for i_th = 1:Threads.nthreads()] # Container for the computation of the amplification
############################################################
# Computes the "dielectric" matrix
# E = I - Psi * M
# ATTENTION, we assume that TAB_PSI has been computed
# @IMPROVE -- The computation of TAB_E could be streamlined a bit
############################################################
function TAB_E!(k::Int64,omega::Complex{Float64},
                tab_E::Matrix{Complex{Float64}})
    #####
    fill!(tab_E,0.0 + 0.0*im) # Initialising with zeros
    #####
    # Adding the identity matrix
    for i=1:NB_INT
        tab_E[i,i] = 1.0 + 0.0*im
    end
    #####
    # Substracting the matrix product PSI * M
    for j=1:NB_INT
        #####
        J = TAB_J[j]
        #####
        val_M = get_M(k,omega,J) # Computing the response matrix
        #####
        for i=1:NB_INT
            tab_E[i,j] -= TAB_PSI[i,j] * val_M # ATTENTION, to the minus sign
        end
        #####
    end
end
############################################################
# Computes the determinant of the dielectric matrix.
# We have found a mode if this quantity goes to zero
# ATTENTION, we compute this quantity for the harmonic k=K
############################################################
function get_detE(k::Int64,omega::Complex{Float64},
                  tab_E::Matrix{Complex{Float64}})
    #####
    TAB_E!(k,omega,tab_E) # Computing the dielectric matrix
    #####
    return det(tab_E) # Returning the determinant
    ##### 
end