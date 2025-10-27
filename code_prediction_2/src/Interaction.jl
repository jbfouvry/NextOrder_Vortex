############################################################
# Computes the radii ra = r_alpha, rb = r_beta
# appearing in the softened pairwise interaction
############################################################
function get_rarb(r::TYPE,rp::TYPE)
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
function UPotk(k::Int64,J::TYPE,Jp::TYPE)
    r, rp = sqrt(2 * J), sqrt(2 * Jp) # Associated radii
    ra, rb = get_rarb(r,rp) # Computing ra, rb
    return (G/(4 * PI)) * (1 // (abs(k))) * (rb/ra)^(abs(k))
end
############################################################
# Gradient of the interaction potential wrt the first
# and second variable
# @IMPROVE -- This expression is not very numerically stable
############################################################
function UPotkP1(k::Int64,J::TYPE,Jp::TYPE)
    r, rp = sqrt(2 * J), sqrt(2 * Jp) # Associated radii
    ra, rb = get_rarb(r,rp) # Computing ra, rb
    drdJ = 1 / r # Gradient dr/dJ
    dUdr = (G / (4 * PI))*
           ((rp^2+EPS^2-r^2)/(r*sqrt(((r+rp)^2+EPS^2)*((r-rp)^2+EPS^2))))*
           (rb/ra)^(abs(k)) # Computing the gradient dU[r,rp]/dr
    #####
    dUdJ = dUdr * drdJ # Gradient dU/dJ
    #####
    return dUdJ # Output
end
#####
function UPotkP2(k::Int64,J::TYPE,Jp::TYPE)
    return UPotkP1(k,Jp,J) # Simply using the symmetry
end
