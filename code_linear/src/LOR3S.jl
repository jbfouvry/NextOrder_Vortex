############################################################
# Mean potential and orbital frequency for the LOR3S profile
# WITHOUT softening
############################################################
function Omega_LOR3S(J::Float64)
    if J <= J_LOR3S # The frequency is zero inside
        return 0.0
    else
        return -0.25*(G*(J-J_LOR3S)*GAMMA_TOT*(J-J_LOR3S+2*SIGMA_LOR3S))/
               (J*pi*(J-J_LOR3S+SIGMA_LOR3S)^(2))
    end
end