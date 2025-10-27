#  Global constants of the problem
const G = 1.0 # Gravitational constant
const GAMMA_TOT = 1.0 # Total circulation
const GAMMA_PART = Q * GAMMA_TOT / NPART # Individual mass of the particles
#####
const DIST_SELF = 10^(-2) # Target error for the convergence of the iterative algorithm
const MAX_ITER_SELF = 10^(6)