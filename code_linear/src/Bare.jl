const TAB_PSI = zeros(Complex{Float64},NB_INT,NB_INT)  # Bare interaction. ATTENTION, this is a complex quantity
############################################################
# Initialising the bare matrix
# @IMPROVE -- In practice, the matrix is symmetric
#             so we could compute only half of the coefficients
############################################################
function TAB_PSI!(k::Int64)
    for j=1:NB_INT
        #####
        J_j = TAB_J[j]
        #####
        for i=1:NB_INT
            #####
            J_i = TAB_J[i]
            #####
            TAB_PSI[i,j] = UPotk(k,J_i,J_j) + 0.0*im # ATTENTION, it is a complex quantity
            #####
        end
    end
end
############################################################
TAB_PSI!(K) # Computing once the bare matrix for the harmonics k=K