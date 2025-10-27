############################################################
# Computes the radii ra = r_alpha, rb = r_beta
# appearing in the softened pairwise interaction
############################################################
function get_rarb(r::Float64,rp::Float64)
    #####
    ra = sqrt((1 // 2) * (r^2 + rp^2 + EPS^2 +
                   sqrt(((r + rp)^2 + EPS^2) * ((r - rp)^2 + EPS^2)))) # Computing r_alpha
    rb = (r * rp) / ra # Computing r_beta
    #####
    return ra, rb
end
############################################################
# Expression of the coupling coefficients in Fourier space
# ATTENTION, we assume |k| >= 1
############################################################
function UPotk(k::Int64,J::Float64,Jp::Float64)
    r, rp = sqrt(2 * J), sqrt(2 * Jp) # Associated radii
    ra, rb = get_rarb(r,rp) # Computing ra, rb
    return (G/(4 * pi)) * (1 // (abs(k))) * (rb/ra)^(abs(k))
end
############################################################
# Definition of the total coupling coefficients,
# summed up to kmax,
# that appear in the Landau equation
############################################################
function UtotSQ(kmax::Int64,J::Float64,Jp::Float64)
    #####
    res = 0.0
    #####
    for k=1:kmax # ATTENTION, we sum only over the k > 0
        res += k * UPotk(k,J,Jp)^(2) # ATTENTION, not yet multiplied by 2
    end
    #####
    res *= 2 # Prefactor
    #####
    return res
end