############################################################
# Definition of the time kernel used for the regularisation
# See Eq. (24) in Bar-Or & Fouvry (2018)
############################################################
function kernel_REG_COS(T::Float64,omega::Float64)
    #####
    x = T * omega
    #####
    # Close to zero, we perform a limited development
    if (abs(x) <= 10^(-3))
        x_SQ = x * x
        res = x_SQ * (-1/(40320 * pi))
        res += 1/(720 * pi)
        res *= x_SQ
        res += (-1/(24 * pi))
        res *= x_SQ
        res += (1/(2 * pi))
        res *= T
        return res
    else
        return (1 - cos(T * omega)) / (pi * omega^2 * T) # Generic expression
    end
end
############################################################
# Regularisation function using a Lorentzian
############################################################
function kernel_REG_LRZ(T::Float64,omega::Float64)
    #####
    return (T / pi) / (1 + (T * omega)^(2))
    #####
end
############################################################
# Wrapped function that computes the rescaled flux
# for a given action, while computing the coupling coefficient
# up to a certain kmax
# ATTENTION, we regularise the DIRAC *after* performing the sum over k
############################################################
const DELTA_J = (2 * SIGMA_LOC) / NB_INT # Step distance to perform the integral
#####
function FluxRESCALED_REG_AFTR(kmax::Int64,T::Float64,J::Float64,
                               kernel_REG::F)  where {F <: Function}
    #####
    res = 0.0 # Initialising the result
    #####
    omega = Omega(J) # Frequency
    #####
    for i = 1:NB_INT
        #####
        Jp = J_LOC - SIGMA_LOC + (i - 0.5) * DELTA_J # Current action
        #####
        omegap = Omega(Jp) # Current frequency
        #####
        res += UtotSQ(kmax,J,Jp) * (DF(Jp) * DFP(J) - DF(J) * DFP(Jp)) *
               kernel_REG(T, omega - Omega(Jp)) # ATTENTION, not yet multiplied by the prefactor
        #####
    end
    #####
    res *= DELTA_J * 2 * pi^(2) # ATTENTION, need to multiply by DeltaJ
    #####
    return res
end
############################################################
# When we regularise each k-scale independently,
# we may provide a different regularising time for each k
############################################################
# No scaling wrt k
############################################################
function get_Tk_0(T::Float64,k::Int64)
    return T
end
############################################################
# Scaling like 1/|k|
############################################################
function get_Tk_1(T::Float64,k::Int64)
    return T / abs(k)
end
############################################################
# Scaling like 1/|k|^2
############################################################
function get_Tk_2(T::Float64, k::Int64)
    return T / abs(k)^2
end
############################################################
# Same as above
# ATTENTION, we regularise the DIRAC *before* performing the sum over k
# Here, we provide the regularising kernel in time
# as well as the scaling function for the non-linear time, as a function of k
############################################################
function FluxRESCALED_REG_BFRE(kmax::Int64,T::Float64,J::Float64,
                               kernel_REG::F1, get_Tk::F2)  where {F1 <: Function, F2 <: Function}
    #####
    res = 0.0 # Initialising the result
    #####
    omega = Omega(J) # Frequency
    #####
    for i = 1:NB_INT
        #####
        Jp = J_LOC - SIGMA_LOC + (i - 0.5) * DELTA_J # Current action
        #####
        omegap = Omega(Jp) # Current frequency
        #####
        diff_omega = omega - Omega(Jp) # Difference in frequency
        #####
        cross_DF = (DF(Jp) * DFP(J) - DF(J) * DFP(Jp)) # Cross term of the DF
        #####
        res_loc = 0.0 # Local counter for the sum over k
        #####
        # By symmetry, k and -k contribute the same
        # So, we will multiply by a factor 2, at the end
        for k=1:KMAX
            #####
            Tk = get_Tk(T,k) # Non-linear for the current scale
            #####
            res_loc += k^(2) * UPotk(k,J,Jp)^(2) * kernel_REG(Tk, k * diff_omega)
        end
        #####
        res += res_loc * cross_DF # Multiplying by the cross DF
    end
    #####
    # Multiplying by the prefactors:
    # + DELTA_J from the midpoint rule
    # + 2 from k <-> -k
    # + 2 * pi^(2) from the prefactor of the kinetic equation
    #####
    res *= DELTA_J * 2 * 2 * pi^(2) # ATTENTION, need to multiply by DeltaJ
    #####
    return res
end
