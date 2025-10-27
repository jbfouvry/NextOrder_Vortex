##################################################
# Initial sampling according to the LOC DF
# of the form F(J) \propto ((J-J_LOC)/SIGMA_LOC + 1)^2 * 
#                          ((J-J_LOC)/SIGMA_LOC - 1)^2 * 
#                          Theta[|J-J_LOC| < SIGMA_LOC]
# with the action J=(1/2)r^2
##################################################
function init_1part_LOC()
    #####
    theta = 2.0*pi*rand() # Random angle in [0,2pi]
    #####
    # We use a rejection sampling
    b = true # Boolean to check the rejection
    J = NaN # Initialising the action
    while b
        #####
        # Drawing a tentative action
        J_try = J_LOC + SIGMA_LOC * (2 * rand() - 1) # Drawing uniformly within the interval [J_LOC - SIGMA_LOC ; J_LOC + SIGMA_LOC]
        #####
        u = rand() # Random number in [0,1]
        #####
        ratio = (((J_try-J_LOC)/SIGMA_LOC) + 1)^(2) * (((J_try-J_LOC)/SIGMA_LOC) - 1)^(2) # Ratio of the PDFs
        #####
        if (u <= ratio)
            b = false
            J = J_try
        end
    end
    #####
    r = sqrt(2.0*J) # Radius of the particle -- J=r^2/2
    #####
    s, c = sincos(theta)
    x, y = r*s, r*c # Position of the particle. ATTENTION, to the convention for the angle. Though, this plays no role for the sampling here
    #####
    return x, y # Output
end
############################################################
# Mean potential for the LOC profile
# WITH softening
# See the MMA notebook "Check.nb"
# ATTENTION, this is expression is NOT numerically stable
# ATTENTION, this does not contain the active fraction Q
############################################################
function U_LOC_soft(J::Float64)
    #####
    v1 = 4*J^(2)+(2*J_LOC+EPS^(2)-2*SIGMA_LOC)^(2)+4*J*(-2*J_LOC+EPS^(2)+2*SIGMA_LOC)
    #####
    v2 = 4*J^(2)+4*J*(-2*J_LOC+EPS^(2)-2*SIGMA_LOC)+(2*J_LOC+EPS^(2)+2*SIGMA_LOC)^(2)
    #####
    res = (G*GAMMATOT*(-3*(sqrt(v1)-sqrt(v2))*(96*J^(4)-3848*J^(3)*EPS^(2)+7524*J^(2)*
    EPS^(4)-2562*J*EPS^(6)+131*EPS^(8))+6*(sqrt(v1)+sqrt(v2))*(48*J^(3)-796*J^(2)*EPS^(2)+
    708*J*EPS^(4)-71*EPS^(6))*SIGMA_LOC+4*(sqrt(v1)-sqrt(v2))*(328*J^(2)-3326*J*EPS^(2)+
    877*EPS^(4))*SIGMA_LOC^(2)-8*(sqrt(v1)+sqrt(v2))*(164*J-337*EPS^(2))*SIGMA_LOC^(3)-
    5888*(sqrt(v1)-sqrt(v2))*SIGMA_LOC^(4)+24*J_LOC^(4)*(-137*(sqrt(v1)-sqrt(v2))+240*
    SIGMA_LOC)+12*J_LOC^(3)*(326*J*(sqrt(v1)-sqrt(v2))-933*sqrt(v1)*EPS^(2)+933*sqrt(v2)*
    EPS^(2)-154*sqrt(v1)*SIGMA_LOC-154*sqrt(v2)*SIGMA_LOC)+2*J_LOC^(2)*(-3*(sqrt(v1)-
    sqrt(v2))*(548*J^(2)-5716*J*EPS^(2)+1507*EPS^(4))+24*(sqrt(v1)+sqrt(v2))*(43*J-94*
    EPS^(2))*SIGMA_LOC+3836*(sqrt(v1)-sqrt(v2))*SIGMA_LOC^(2)-8640*SIGMA_LOC^(3))+J_LOC*
    (1512*J^(3)*(sqrt(v1)-sqrt(v2))-3099*sqrt(v1)*EPS^(6)+3099*sqrt(v2)*EPS^(6)-2466*
    sqrt(v1)*EPS^(4)*SIGMA_LOC-2466*sqrt(v2)*EPS^(4)*SIGMA_LOC+13012*sqrt(v1)*EPS^(2)*
    SIGMA_LOC^(2)-13012*sqrt(v2)*EPS^(2)*SIGMA_LOC^(2)+3352*sqrt(v1)*SIGMA_LOC^(3)+3352*
    sqrt(v2)*SIGMA_LOC^(3)-36*J^(2)*(959*(sqrt(v1)-sqrt(v2))*EPS^(2)+34*(sqrt(v1)+sqrt(v2))*
    SIGMA_LOC)+J*(31602*(sqrt(v1)-sqrt(v2))*EPS^(4)+9408*(sqrt(v1)+sqrt(v2))*EPS^(2)*
    SIGMA_LOC-4664*(sqrt(v1)-sqrt(v2))*SIGMA_LOC^(2)))+1024*SIGMA_LOC^(5)*(23+30*log(2))+
    60*((3*(16*J_LOC^(5)+80*(J-J_LOC)^(4)*EPS^(2)+80*(J-J_LOC)^(2)*(-4*J+J_LOC)*EPS^(4)+40*
    (6*J^(2)-6*J*J_LOC+J_LOC^(2))*EPS^(6)+10*(-4*J+J_LOC)*EPS^(8)+EPS^(10))-40*(4*J_LOC^(3)+
    12*(J-J_LOC)^(2)*EPS^(2)+6*(-2*J+J_LOC)*EPS^(4)+EPS^(6))*SIGMA_LOC^(2)+240*(J_LOC+EPS^(2))*
    SIGMA_LOC^(4)+128*SIGMA_LOC^(5))*log(-2*J+2*J_LOC+sqrt(v1)+EPS^(2)-2*SIGMA_LOC)-3*(16*
    J_LOC^(5)+80*(J-J_LOC)^(4)*EPS^(2)+80*(J-J_LOC)^(2)*(-4*J+J_LOC)*EPS^(4)+40*(6*J^(2)-6*
    J*J_LOC+J_LOC^(2))*EPS^(6)+10*(-4*J+J_LOC)*EPS^(8)+EPS^(10))*log(-2*J+2*J_LOC+sqrt(v2)+
    EPS^(2)+2*SIGMA_LOC)+8*(-64*SIGMA_LOC^(5)*log(2*J+2*J_LOC+sqrt(v1)+EPS^(2)-2*SIGMA_LOC)+
    SIGMA_LOC^(2)*(5*(4*J_LOC^(3)+12*(J-J_LOC)^(2)*EPS^(2)+6*(-2*J+J_LOC)*EPS^(4)+EPS^(6))-
    30*(J_LOC+EPS^(2))*SIGMA_LOC^(2)-16*SIGMA_LOC^(3))*log(-2*J+2*J_LOC+sqrt(v2)+EPS^(2)+2*
    SIGMA_LOC)+2*(J_LOC+SIGMA_LOC)^(3)*(3*J_LOC^(2)-9*J_LOC*SIGMA_LOC+8*SIGMA_LOC^(2))*(-
    log(4*J^(2)+2*J*(-2*J_LOC+sqrt(v2)+2*EPS^(2)-2*SIGMA_LOC)+EPS^(2)*(2*J_LOC+sqrt(v2)+
    EPS^(2)+2*SIGMA_LOC))+log(4*J^(2)+EPS^(2)*(2*J_LOC+sqrt(v1)+EPS^(2)-2*SIGMA_LOC)+2*J*
    (-2*J_LOC+sqrt(v1)+2*(EPS^(2)+SIGMA_LOC))))))))/(122880*pi*SIGMA_LOC^(5))
    #####
    return res
end
############################################################
# Mean orbital frequency for the LOC profile
# ATTENTION, this expression is NOT numerically stable
# ATTENTION, this does not contain the active fraction Q
############################################################
@fastmath function Omega_LOC_soft(J::Float64)
    #####
    v1 = 4*J^(2)+(2*J_LOC+EPS^(2)-2*SIGMA_LOC)^(2)+4*J*(-2*J_LOC+EPS^(2)+2*SIGMA_LOC)
    #####
    v2 = 4*J^(2)+4*J*(-2*J_LOC+EPS^(2)-2*SIGMA_LOC)+(2*J_LOC+EPS^(2)+2*SIGMA_LOC)^(2)
    #####
    res = (G*GAMMATOT*(-((sqrt(v1)-sqrt(v2))*(48*(J-J_LOC)^(4)-16*(89*J-6*J_LOC)*(J-J_LOC)^(2)*
    EPS^(2)+8*(239*J^(2)-184*J*J_LOC+9*J_LOC^(2))*EPS^(4)+4*(-89*J+6*J_LOC)*EPS^(6)+3*EPS^(8)))+
    2*(sqrt(v1)+sqrt(v2))*(2*J-2*J_LOC-EPS^(2))*(12*(J-J_LOC)^(2)+4*(-32*J+3*J_LOC)*EPS^(2)+3*
    EPS^(4))*SIGMA_LOC+4*(sqrt(v1)-sqrt(v2))*(28*(J-J_LOC)^(2)+4*(-39*J+7*J_LOC)*EPS^(2)+7*
    EPS^(4))*SIGMA_LOC^(2)-56*(sqrt(v1)+sqrt(v2))*(2*J-2*J_LOC-EPS^(2))*SIGMA_LOC^(3)+128*(-
    sqrt(v1)+sqrt(v2))*SIGMA_LOC^(4)-512*SIGMA_LOC^(5)+240*J*EPS^(2)*(2*J-2*J_LOC-EPS^(2))*(4*
    (J-J_LOC)^(2)+2*(-5*J+2*J_LOC)*EPS^(2)+EPS^(4)-4*SIGMA_LOC^(2))*(log(2*J-2*J_LOC+sqrt(v2)-
    EPS^(2)-2*SIGMA_LOC)-log(2*J-2*J_LOC+sqrt(v1)-EPS^(2)+2*SIGMA_LOC))))/(4096*J*pi*SIGMA_LOC^(5))
    #####
    return res
end