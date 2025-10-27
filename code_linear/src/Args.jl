##################################################
# Parsing of the command-line arguments
##################################################
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--eps"
    help = "Softening length"
    arg_type = Float64
    default = 0.01
    "--q"
    help = "Active fraction"
    arg_type = Float64
    default = 0.1
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
    "--k"
    help = "Considered harmonics"
    arg_type = Int64
    default = 1
    "--T"
    help = "Time for the Hilbert kernel -- T=0.0 is no regularisation"
    arg_type = Float64
    default = 10.0
    "--nb_int"
    help = "Number of integration points"
    arg_type = Int64
    default = 1000
end
##################################################
parsed_args = parse_args(tabargs)
##################################################
# General parameters
##################################################
const EPS = parsed_args["eps"] # Softening length
const Q = parsed_args["q"] # Active fraction
#####
const BACKGROUND = parsed_args["background"] # Background considered
#####
if ((BACKGROUND != "LOR3") && (BACKGROUND != "LOR3S"))
    error("Supported BACKGROUND: true/false")
end
#####
const SIGMA_LOR3 = parsed_args["sigma_LOR3"] # Spread of the background
const J_LOC = parsed_args["J_LOC"] # Central action for the LOC DF
const SIGMA_LOC = parsed_args["sigma_LOC"] # Width for the LOC DF
const SIGMA_LOR3S = parsed_args["sigma_LOR3S"] # Spread of the LOR3S background
const J_LOR3S = parsed_args["J_LOR3S"] # Offset of the LOR3S background
const K = parsed_args["k"] # Considered harmonics
const T = parsed_args["T"] # Width of the Hilbert kernel
const NB_INT = parsed_args["nb_int"] # Number of integration points
##################################################
# Sanity checks
##################################################
if (J_LOR3S >= SIGMA_LOR3S)
    error("J_LOR3S too large compared to SIGMA_LOR3S")
end
#####
if (J_LOC < SIGMA_LOC)
    error("J_LOC should be larger than SIGMA_LOC")
end