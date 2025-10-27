############################################################
# Definition of the Hilbert kernel
# 1/z_{-}, with regularisation
# @IMPROVE -- We should perform a limited development for t * omega -> 0
############################################################
function get_Hilbert_regular(z::Complex{Float64},t::Float64)
    return (1.0 - exp(- im * z * t)) / (z)
end
############################################################
# Hilbert kernel without regularisation
# ATTENTION, this is to be used only in the upper half of the complex plane
############################################################
function get_Hilbert_upper(z::Complex{Float64},t::Float64)
    return 1.0 / z
end
############################################################
# If the provided T is 0.0 or negative, we pick the "upper" Hilbert kernel
if (T <= 0)
    const get_Hilbert = get_Hilbert_upper # Without regularisation
else
    const get_Hilbert = get_Hilbert_regular
end

