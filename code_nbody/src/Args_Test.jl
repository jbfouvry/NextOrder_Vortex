##################################################
# Parsing of the command-line arguments
##################################################
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--Npart"
    help = "Total number of particles"
    arg_type = Int64
    default = 16
    "--Ntest"
    help = "Number of zero-mass test particles"
    arg_type = Int64
    default = 100
    "--Jtest"
    help = "Initial action of the test particles"
    arg_type = Float64
    default = 1.0
    "--q"
    help = "Active fraction"
    arg_type = Float64
    default = 1.0
    "--eps"
    help = "Softening length"
    arg_type = Float64
    default = 0.05
    "--recenter"
    help = "Whether to recenter the initial sampling: true/false"
    arg_type = String
    default = "false"
    "--background"
    help = "Background considered: LOR3/LOR3S"
    arg_type = String
    default = "LOR3"
    "--sigma_LOR3"
    help = "Spread of the LOR3 background"
    arg_type = Float64
    default = 1.0
    "--sigma_LOR3S"
    help = "Spread of the LOR3S background"
    arg_type = Float64
    default = 1.0
    "--J_LOR3S"
    help = "Offset for the LOR3S profile"
    arg_type = Float64
    default = 0.5
    "--J_LOC"
    help = "Central action for the LOC DF"
    arg_type = Float64
    default = 1.0
    "--sigma_LOC"
    help = "Width for the LOC DF"
    arg_type = Float64
    default = 0.2
    "--scheme"
    help = "Integration scheme: RK1/RK2/RK4/RK6/RK9/RK14"
    arg_type = String
    default = "RK4"
    "--dt"
    help = "Integration timestep"
    arg_type = Float64
    default = 0.1
    "--Nsteps"
    help = "Number of steps dt used in integrate_Nsteps!()"
    arg_type = Int64
    default = 1000
    "--Ndumps"
    help = "Number of dumps"
    arg_type = Int64
    default = 100
    "--seed"
    help = "Seed for the random generator"
    arg_type = Int64
    default = 1
end
##################################################
parsed_args = parse_args(tabargs)
##################################################
# General parameters
##################################################
const NPART = parsed_args["Npart"]  # Total number of particles
const NTEST = parsed_args["Ntest"]  # Number of zero-mass test particles
const NMASS = NPART - NTEST # Number of massive particles
#####
if (NMASS <= 0)
    error("NMASS should be positive")
end
#####
const JTEST = parsed_args["Jtest"] # Initial action of the test particles
const Q = parsed_args["q"] # Self-gravitating fraction
const EPS = parsed_args["eps"] # Softening length
#####
const RECENTER = parsed_args["recenter"] # Whether to recenter the initial sampling
#####
if ((RECENTER != "true") && (RECENTER != "false"))
    error("Supported RECENTER: true/false")
end
#####
const BACKGROUND = parsed_args["background"] # Background considered
#####
if ((BACKGROUND != "LOR3") && (BACKGROUND != "LOR3S"))
    error("Supported BACKGROUND: true/false")
end
#####
const SIGMA_LOR3 = parsed_args["sigma_LOR3"] # Spread of the LOR3 background
const SIGMA_LOR3S = parsed_args["sigma_LOR3S"] # Spread of the LOR3S background
const J_LOR3S = parsed_args["J_LOR3S"] # Offset of the LOR3S background
const J_LOC = parsed_args["J_LOC"] # Central action for the LOC DF
const SIGMA_LOC = parsed_args["sigma_LOC"] # Width for the LOC DF
#####
if (J_LOR3S >= SIGMA_LOR3S)
    error("J_LOR3S too large compared to SIGMA_LOR3S")
end
#####
if (J_LOC < SIGMA_LOC)
    error("SIGMA_LOC too large compared to J_LOC")
end
#####
const SCHEME = parsed_args["scheme"] # Integration scheme
#####
if ((SCHEME != "RK1") &&
    (SCHEME != "RK2") &&
    (SCHEME != "RK4") &&
    (SCHEME != "RK6") &&
    (SCHEME != "RK9") &&
    (SCHEME != "RK14"))
    error("Supported SCHEME: RK1/RK2/RK4/RK6/RK9/RK14")
end
#####
const DT     = parsed_args["dt"]     # Integration timestep
const NSTEPS = parsed_args["Nsteps"] # Number of integration timesteps between every dump
const NDUMPS = parsed_args["Ndumps"] # Number of dumps of data
const SEED   = parsed_args["seed"]   # Sampling seed
##################################################
Random.seed!(SEED) # Setting the random seed
##################################################
