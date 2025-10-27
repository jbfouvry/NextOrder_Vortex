##################################################
# Parsing of the command-line arguments
##################################################
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--Npart"
    help = "Number of particles"
    arg_type = Int64
    default = 1000
    "--eps"
    help = "Softening length"
    arg_type = Float64
    default = 0.01
    "--q"
    help = "Active fraction"
    arg_type = Float64
    default = 0.0001
    "--sigma_LOR3S"
    help = "Spread of the LOR3S background"
    arg_type = Float64
    default = 1.0
    "--J_LOR3S"
    help = "Location of the LOR3S background"
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
    "--kmax"
    help = "Fundamental resonance considered"
    arg_type = Int64
    default = 100
    "--nb_int"
    help = "Number of integration points for the regularisation in Flux_REG"
    arg_type = Int64
    default = 1000
    "--Treg_over_Tdyn"
    help = "Regularisation time in Flux_REG over TDYN"
    arg_type = Float64
    default = 10.0
end
##################################################
parsed_args = parse_args(tabargs)
##################################################
# General parameters
##################################################
const NPART = parsed_args["Npart"] # Total number of particles
const EPS = parsed_args["eps"] # Softening length
const Q = parsed_args["q"] # Active fraction
const SIGMA_LOR3S = parsed_args["sigma_LOR3S"] # Spread of the background
const J_LOR3S = parsed_args["J_LOR3S"] # Position of the background
const J_LOC = parsed_args["J_LOC"] # Central action for the LOC DF
const SIGMA_LOC = parsed_args["sigma_LOC"] # Width for the LOC DF
const KMAX = parsed_args["kmax"] # Maximum range for the fundamental resonances
const NB_INT = parsed_args["nb_int"] # Number of integration points in Flux_REG
const TREG_OVER_TDYN = parsed_args["Treg_over_Tdyn"] # Regularisation time in Flux_REG over Tdyn
#####
if (J_LOR3S >= SIGMA_LOR3S)
    error("J_LOR3S too large compared to SIGMA_LOR3S")
end
#####
if (J_LOC < SIGMA_LOC)
    error("SIGMA_LOC too large compared to J_LOC")
end
#####