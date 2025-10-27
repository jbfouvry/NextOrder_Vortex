############################################################
# Expression of the coupling coefficients that appear
# in the kinetic equation
############################################################
# First coupling coefficient U1
############################################################
function get_U1(k1::Int64,k2::Int64,
                J::TYPE,J1::TYPE,J2::TYPE)
    #####
    res = k2*(k1+k2)*
          (UPotk(k1+k2,J,J2)*UPotkP2(k1,J1,J2) -
           UPotk(k2,J,J2)   *UPotkP1(k1,J,J1)) +
           k1*(k1+k2)*
          (UPotk(k1,J,J1)   *UPotkP1(k2,J,J2) -
           UPotk(k1+k2,J,J1)*UPotkP1(k2,J1,J2)) -
           k1*k2*
          (UPotk(k2,J1,J2)*UPotkP2(k1+k2,J,J1) -
           UPotk(k1,J1,J2)*UPotkP2(k1+k2,J,J2))
    #####
    return res
end
############################################################
# For the bar contributions
# This simply amounts to the substitution J1 <-> J2
############################################################
function get_U1bar(k1::Int64,k2::Int64,
                   J::TYPE,J1::TYPE,J2::TYPE)
    #####
    return get_U1(k1,k2,J,J2,J1) # Attention to the order of the arguments
end
############################################################
# Second coupling coefficient U2
############################################################
function get_U2(k1::Int64,k2::Int64,
                J::TYPE,J1::TYPE,J2::TYPE)
    #####
    res = (k1+k2)*OmegaP(J) *UPotk(k1,J,J1)   *UPotk(k2,J,J2) -
           k1    *OmegaP(J1)*UPotk(k1+k2,J,J1)*UPotk(k2,J1,J2) -
           k2    *OmegaP(J2)*UPotk(k1,J1,J2)  *UPotk(k1+k2,J,J2)
    #####
    return res
end
############################################################
# For the bar contributions
# This simply amounts to the substitution J1 <-> J2 in get_U2
############################################################
function get_U2bar(k1::Int64,k2::Int64,
                   J::TYPE,J1::TYPE,J2::TYPE)
    #####
    return get_U2(k1,k2,J,J2,J1) # Attention to the order of the arguments
end
############################################################
# Full expression of the integrand C
# ATTENTION, does not contain yet the resonance condition
# DIRAC[(k1+k2)*Omega[J]-k1*Omega[J1]-k2*Omega[J2]]
############################################################
function get_C(k1::Int64,k2::Int64,
               J::TYPE,J1::TYPE,J2::TYPE)
    #####
    # Computing the crossed term of derivatives
    # ATTENTION, for LOC, we do not use the ratio F'[J]/F[J]
    # crossedDF = DF(J) * DF(J1) * DF(J2) * (
    #             (k1+k2) * DFPoverDF(J)  -
    #              k1     * DFPoverDF(J1) -
    #              k2     * DFPoverDF(J2))
    crossedDF = (k1+k2) * DFP(J) *  DF(J1) *  DF(J2) -
                 k1     *  DF(J) * DFP(J1) *  DF(J2) -
                 k2     *  DF(J) *  DF(J1) * DFP(J2)
    #####
    res = (1 // (k1^(2) * (k1+k2))) *
          (get_U1(k1,k2,J,J1,J2) * diffOmega(J,J1) +
           k2 * get_U2(k1,k2,J,J1,J2))^(2) *
           crossedDF
    #####
    return res
end
#####
# For the bar contributions
# The integrand has the shape
# ATTENTION, does not contain yet the resonance condition
# DIRAC[(k1+k2)*Omega[J]-k1*Omega[J2]-k2*Omega[J1]]
# !! ATTENTION !! the 2nd and 3rd terms are
# - k2 F1'/F1 - k1 F2'/F2
#####
function get_Cbar(k1::Int64,k2::Int64,
                  J::TYPE,J1::TYPE,J2::TYPE)
    #####
    # Computing the crossed term of derivatives
    # ATTENTION, for LOC, we do not use the ratio F'[J]/F[J]
    # crossedDF = DF(J) * DF(J1) * DF(J2) * (
    #             (k1+k2) * DFPoverDF(J)  -
    #              k2     * DFPoverDF(J1) - # ATTENTION, this is "k2"
    #              k1     * DFPoverDF(J2))  # ATTENTION, this is "k1"
    # ATTENTION, to the "k1" and "k2" and their positions
    crossedDF = (k1+k2) * DFP(J) *  DF(J1) *  DF(J2) -
                 k2     *  DF(J) * DFP(J1) *  DF(J2) - # ATTENTION, this is "k2"
                 k1     *  DF(J) *  DF(J1) * DFP(J2)   # ATTENTION, this is "k1"
    #####
    res = (1 // (k2^(2) * (k1+k2))) *
          (get_U1bar(k1,k2,J,J1,J2) * diffOmega(J,J1)-
           k1 * get_U2bar(k1,k2,J,J1,J2))^(2) *
           crossedDF
    #####
    return res
end
############################################################
# Performing the integral over J2
############################################################
# The integrand has the shape
# INT[get_C*DIRAC[(k1+k2)*Omega[J]-k1*Omega[J1]-k2*Omega[J2]],{J2}]
#####
function calc_C(k1::Int64,k2::Int64,
                J::TYPE,J1::TYPE)
    #####
    omg = ((k1+k2) * Omega(J) - k1 * Omega(J1)) / (k2) # Resonance frequency
    #####
    J2_res = Jres(omg) # Determining the resonance frequency
    #####
    if (J2_res != J2_res) # There are no resonance. ATTENTION, this is a lazy way to test for NaN
        return TYPE(0)
    end
    #####
    # Computing the contribution. ATTENTION, not to forget the denominators
    return get_C(k1,k2,J,J1,J2_res) / (abs(k2)*abs(OmegaP(J2_res)))
end
#####
# The integrand has the shape
# INT[get_Cbar*DIRAC[(k1+k2)*Omega[J]-k1*Omega[J2]-k2*Omega[J1]],{J2}]
#####
function calc_Cbar(k1::Int64,k2::Int64,
                   J::TYPE,J1::TYPE)
    #####
    omg = ((k1+k2) * Omega(J) - k2 * Omega(J1)) / (k1) # Resonance frequency
    #####
    J2_res = Jres(omg) # Determining the resonance frequency
    #####
    if (J2_res != J2_res) # There are no resonance. ATTENTION, this is a lazy way to test for NaN
        return TYPE(0)
    end
    #####
    # Computing the contribution. ATTENTION, not to forget the denominators
    return get_Cbar(k1,k2,J,J1,J2_res)/(abs(k1)*abs(OmegaP(J2_res)))
end
############################################################
# The diffusion equation is generically written as
# dF[J]/dt = 2*pi^3*gamma^2*d/dJ[SUM[FluxFund[{k,kp},J],{k,kp}]]
# with
# FluxFund[{k,kp},J] = PRPA[INT[1/((Omega[J]-Omega[J1])^(4))]*
#                      *WellPosedFund[{k,kp},J,J1]]
# where WellPosedFund[{k,kp},J,J1] is the well-posed integrand
# devised by using appropriately tailored averaging weights
############################################################
# Weight used to obtain a well-posed integrand
############################################################
function get_gamma_weight(k::Int64,kp::Int64)
    return TYPE(((k - kp)^(2)) // (k^(2) + kp^(2)))
end
############################################################
# Well-posed expression of the integrand for a given
# fundamental resonance
############################################################
function WellPosedFund(k::Int64,kp::Int64,
                       J::TYPE,J1::TYPE)
    #####
    gamma_loc = get_gamma_weight(k,kp) # Weigth used for the average
    #####
    res = TYPE(0) # Initialising the result
    #####
    # Contribution from (k1,k2)=(k,kp)
    res += calc_C(k,kp,J,J1)
    #####
    # Contribution from (k1,k2)=(kp,k)
    res += calc_C(kp,k,J,J1)
    #####
    # Contribution from (k1,k2)=(k+kp,-k)
    res += calc_C(k+kp,-k,J,J1)
    #####
    # Contribution from (k1,k2)=(k+kp,-kp)
    res += calc_C(k+kp,-kp,J,J1)
    #####
    # Contribution from (k1,k2)=(k,-k-kp)
    res += gamma_loc*calc_C(k,-k-kp,J,J1) + (1 - gamma_loc)*calc_Cbar(k,-k-kp,J,J1)
    #####
    # Contribution from (k1,k2)=(kp,-k-kp)
    res += gamma_loc*calc_C(kp,-k-kp,J,J1) + (1 - gamma_loc)*calc_Cbar(kp,-k-kp,J,J1)
    #####
    return res # Output
end
############################################################
# We have
# dF[J]/dt = SUM[2*pi^3*gamma^2*d/dJ[INT[1/((Omega[J]-Omega[J1])^4) WellPosedFund[k,kp,J,J1],{J1}]],
#            {k,kp}]
# We define the rescaled flux as
# FluxRESCALED[k,kp,J] = SUM[2*pi^3*INT[1/((Omega[J]-Omega[J1])^4) WellPosedFund[k,kp,J,J1],{J1}],
#                        {k,kp}],
# i.e. it does not contain the prefactor gamma^2
# We have dF[J]/dt=gamma^2*d/dJ[SUM[FluxRESCALED[k,kp,J],{k,kp}]]
# #####
# We rewrite the integrand of FluxRESCALED as
# FluxRESCALED[k,kp,J] = INT[IntFluxRESCALED[k,kp,J,J1],{J1}]
# with
# IntFluxRESCALED[k,kp,J,J1] = 2*pi^3*1/((Omega[J]-Omega[J1])^4) WellPosedFund[k,kp,J,J1]
############################################################
function IntFluxRESCALED(k::Int64,kp::Int64,
                         J::TYPE,J1::TYPE)
    #####
    return 2 * PI^(3) * (diffOmega(J,J1))^(-4) * WellPosedFund(k,kp,J,J1)
end
