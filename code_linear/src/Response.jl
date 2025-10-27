############################################################
# Returns the value of the response matrix
# ATTENTION, omega is complex
############################################################
function get_M(k::Int64,omega::Complex{Float64},J::Float64)
    #####
    val = 2 * pi * DELTA_J * (k * DFP(J)) * get_Hilbert(k * Omega(J) - omega,T)
    #####
    return val
end