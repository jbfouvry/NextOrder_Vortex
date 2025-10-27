##################################################
# Parsing of the command-line arguments
##################################################
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--type"
    help = "Floating-point type: Float64/Double64"
    arg_type = String
    default = "Float64"
    "--eps"
    help = "Softening length"
    arg_type = Float64
    default = 0.05
    "--q"
    help = "Active fraction"
    arg_type = Float64
    default = 0.1
    "--sigma_LOR3"
    help = "Spread of the LOR3 background"
    arg_type = Float64
    default = 1.0
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
    default = 50
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
# Parsing the integration type
const TYPE_READ = parsed_args["type"]
if     (TYPE_READ == "Float64")
    const TYPE = Float64
elseif (TYPE_READ == "Double64")
    using DoubleFloats
    const TYPE = Double64
else
    error("Unsupported TYPE")
end
##################################################
# ATTENTION, we parse the input arguments from Float64 to TYPE
# This (lazy) approach has the drawback that, e.g.,
# eps = 0.1 --> 1 + EPS = 1.10000000000000000555111512312578270
# i.e. we lose precision in the definition of the arguments
# But this does not affect the upcoming calculation [I believe].\
# Reading the command-line arguments in Float64
const EPS_READ = parsed_args["eps"] # Softening length
const Q_READ = parsed_args["q"] # Active fraction
const SIGMA_LOR3_READ = parsed_args["sigma_LOR3"] # Spread of the background
const J_LOC_READ = parsed_args["J_LOC"] # Central action for the LOC DF
const SIGMA_LOC_READ = parsed_args["sigma_LOC"] # Width for the LOC DF
#####
# Now, we cast all these arguments to TYPE
const EPS = TYPE(EPS_READ)
const Q = TYPE(Q_READ)
const SIGMA_LOR3 = TYPE(SIGMA_LOR3_READ)
const J_LOC = TYPE(J_LOC_READ)
const SIGMA_LOC =TYPE(SIGMA_LOC_READ)
#####
# These arguments are Int64, so no need to be careful
const KMAX = parsed_args["kmax"] # Maximum range for the fundamental resonances
const NB_INT = parsed_args["nb_int"] # Number of integration points
