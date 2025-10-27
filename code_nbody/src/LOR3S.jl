############################################################
# Mean potential and orbital frequency for the LOR3S profile
# WITHOUT softening
############################################################
# Constant value for the Hamiltonian for J < J_LOR3S
const H0_LOR3S = (G*GAMMATOT*(2*J_LOR3S*SIGMA_LOR3S*log(J_LOR3S)-(J_LOR3S)^(2)*
                  log(2*J_LOR3S)+SIGMA_LOR3S*(SIGMA_LOR3S+J_LOR3S*(-1+log(4))-SIGMA_LOR3S*
                  log(2*SIGMA_LOR3S))))/(4*pi*(J_LOR3S-SIGMA_LOR3S)^(2))
############################################################
function U_LOR3S(J::Float64)
    if J <= J_LOR3S # The Hamiltonian is constant inside
        return H0_LOR3S
    else
        return -0.25*(G*GAMMATOT*((J_LOR3S-SIGMA_LOR3S)*SIGMA_LOR3S^(2)*
               (J+(J_LOR3S-SIGMA_LOR3S)*(-1+log(2)))+(J_LOR3S-SIGMA_LOR3S)^(2)*
               (SIGMA_LOR3S^(2)*log(J)+(J-J_LOR3S)*(J-J_LOR3S+2*SIGMA_LOR3S)*
               log(2*J))+SIGMA_LOR3S^(2)*(J-J_LOR3S+SIGMA_LOR3S)^(2)*
               log((J-J_LOR3S+SIGMA_LOR3S)/J)))/(pi*(J_LOR3S-SIGMA_LOR3S)^(2)*
               (J-J_LOR3S+SIGMA_LOR3S)^(2))
    end
end
############################################################
function Omega_LOR3S(J::Float64)
    if J <= J_LOR3S # The frequency is zero inside
        return 0.0
    else
        return -0.25*(G*(J-J_LOR3S)*GAMMATOT*(J-J_LOR3S+2*SIGMA_LOR3S))/
               (J*pi*(J-J_LOR3S+SIGMA_LOR3S)^(2))
    end
end