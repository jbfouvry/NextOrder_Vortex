# Temporary containers for the integrations
const TAB_X_1 = zeros(Float64,NPART)
const TAB_Y_1 = zeros(Float64,NPART)
#####
const TAB_XDOT_1 = zeros(Float64,NPART)
const TAB_YDOT_1 = zeros(Float64,NPART)
const TAB_XDOT_2 = zeros(Float64,NPART)
const TAB_YDOT_2 = zeros(Float64,NPART)
const TAB_XDOT_3 = zeros(Float64,NPART)
const TAB_YDOT_3 = zeros(Float64,NPART)
const TAB_XDOT_4 = zeros(Float64,NPART)
const TAB_YDOT_4 = zeros(Float64,NPART)
const TAB_XDOT_5 = zeros(Float64,NPART)
const TAB_YDOT_5 = zeros(Float64,NPART)
const TAB_XDOT_6 = zeros(Float64,NPART)
const TAB_YDOT_6 = zeros(Float64,NPART)
const TAB_XDOT_7 = zeros(Float64,NPART)
const TAB_YDOT_7 = zeros(Float64,NPART)
const TAB_XDOT_8 = zeros(Float64,NPART)
const TAB_YDOT_8 = zeros(Float64,NPART)
###############################################
# Pre-computed coefficients for the RK6 integrator
# This follows from Table 5.4. in Hairer
# "Solving Ordinary Differential Equations"
# This is Verner's method of order 6(5)
###############################################
const b1_RK6 = convert(Float64, 3 // 40)
const b3_RK6 = convert(Float64, 875 // 2244)
const b4_RK6 = convert(Float64, 23 // 72)
const b5_RK6 = convert(Float64, 264 // 1955)
const b7_RK6 = convert(Float64, 125 // 11592)
const b8_RK6 = convert(Float64, 43 // 616)

const a10_RK6 = convert(Float64, 1 // 6)

const a20_RK6 = convert(Float64, 4 // 75)
const a21_RK6 = convert(Float64, 16 // 75)

const a30_RK6 = convert(Float64, 5 // 6)
const a31_RK6 = convert(Float64, -8 // 3)
const a32_RK6 = convert(Float64, 5 // 2)

const a40_RK6 = convert(Float64, -165 // 64)
const a41_RK6 = convert(Float64, 55 // 6)
const a42_RK6 = convert(Float64, -425 // 64)
const a43_RK6 = convert(Float64, 85 // 96)

const a50_RK6 = convert(Float64, 12 // 5)
const a51_RK6 = convert(Float64, -8)
const a52_RK6 = convert(Float64, 4015 // 612)
const a53_RK6 = convert(Float64, -11 // 36)
const a54_RK6 = convert(Float64, 88 // 255)

const a60_RK6 = convert(Float64, -8263 // 15000)
const a61_RK6 = convert(Float64, 124 // 75)
const a62_RK6 = convert(Float64, -643 // 680)
const a63_RK6 = convert(Float64, -81 // 250)
const a64_RK6 = convert(Float64, 2484 // 10625)

const a70_RK6 = convert(Float64, 3501 // 1720)
const a71_RK6 = convert(Float64, -300 // 43)
const a72_RK6 = convert(Float64, 297275 // 52632)
const a73_RK6 = convert(Float64, -319 // 2322)
const a74_RK6 = convert(Float64, 24068 // 84065)
const a76_RK6 = convert(Float64, 3850 // 26703)
##################################################
# Integrates for one timestep using RK6
##################################################
function integrate_DT_RK6!()
    ####################
    copy_tab!(TAB_X,TAB_Y,TAB_X_1,TAB_Y_1) # Initial position
    ####################
    # Step 1
    ####################
    TAB_DOT!()
    copy_tab!(TAB_XDOT,TAB_YDOT,TAB_XDOT_1,TAB_YDOT_1)
    ####################
    # Step 2
    ####################
    combine_tab!((a10_RK6,),
                 (TAB_XDOT_1,),
                 (TAB_YDOT_1,),
                  TAB_XDOT,
                  TAB_YDOT)
    drift!(TAB_X_1,TAB_Y_1,
           TAB_X,TAB_Y,
           TAB_XDOT,TAB_YDOT,
           DT)
    TAB_DOT!()
    copy_tab!(TAB_XDOT,TAB_YDOT,TAB_XDOT_2,TAB_YDOT_2)
    ####################
    # Step 3
    ####################
    combine_tab!((a20_RK6,a21_RK6),
                 (TAB_XDOT_1,TAB_XDOT_2),
                 (TAB_YDOT_1,TAB_YDOT_2),
                  TAB_XDOT,
                  TAB_YDOT)
    drift!(TAB_X_1,TAB_Y_1,
           TAB_X,TAB_Y,
           TAB_XDOT,TAB_YDOT,
           DT)
    TAB_DOT!()
    copy_tab!(TAB_XDOT,TAB_YDOT,TAB_XDOT_3,TAB_YDOT_3)
    ####################
    # Step 4
    ####################
    combine_tab!((a30_RK6,a31_RK6,a32_RK6),
                 (TAB_XDOT_1,TAB_XDOT_2,TAB_XDOT_3),
                 (TAB_YDOT_1,TAB_YDOT_2,TAB_YDOT_3),
                  TAB_XDOT,
                  TAB_YDOT)
    drift!(TAB_X_1,TAB_Y_1,
           TAB_X,TAB_Y,
           TAB_XDOT,TAB_YDOT,
           DT)
    TAB_DOT!()
    copy_tab!(TAB_XDOT,TAB_YDOT,TAB_XDOT_4,TAB_YDOT_4)
    ####################
    # Step 5
    ####################
    combine_tab!((a40_RK6,a41_RK6,a42_RK6,a43_RK6),
                 (TAB_XDOT_1,TAB_XDOT_2,TAB_XDOT_3,TAB_XDOT_4),
                 (TAB_YDOT_1,TAB_YDOT_2,TAB_YDOT_3,TAB_YDOT_4),
                  TAB_XDOT,
                  TAB_YDOT)
    drift!(TAB_X_1,TAB_Y_1,
           TAB_X,TAB_Y,
           TAB_XDOT,TAB_YDOT,
           DT)
    TAB_DOT!()
    copy_tab!(TAB_XDOT,TAB_YDOT,TAB_XDOT_5,TAB_YDOT_5)
    ####################
    # Step 6
    ####################
    combine_tab!((a50_RK6,a51_RK6,a52_RK6,a53_RK6,a54_RK6),
                 (TAB_XDOT_1,TAB_XDOT_2,TAB_XDOT_3,TAB_XDOT_4,TAB_XDOT_5),
                 (TAB_YDOT_1,TAB_YDOT_2,TAB_YDOT_3,TAB_YDOT_4,TAB_YDOT_5),
                  TAB_XDOT,
                  TAB_YDOT)
    drift!(TAB_X_1,TAB_Y_1,
           TAB_X,TAB_Y,
           TAB_XDOT,TAB_YDOT,
           DT)
    TAB_DOT!()
    copy_tab!(TAB_XDOT,TAB_YDOT,TAB_XDOT_6,TAB_YDOT_6)
    ####################
    # Step 7
    ####################
    combine_tab!((a60_RK6,a61_RK6,a62_RK6,a63_RK6,a64_RK6),
                 (TAB_XDOT_1,TAB_XDOT_2,TAB_XDOT_3,TAB_XDOT_4,TAB_XDOT_5),
                 (TAB_YDOT_1,TAB_YDOT_2,TAB_YDOT_3,TAB_YDOT_4,TAB_YDOT_5),
                  TAB_XDOT,
                  TAB_YDOT)
    drift!(TAB_X_1,TAB_Y_1,
           TAB_X,TAB_Y,
           TAB_XDOT,TAB_YDOT,
           DT)
    TAB_DOT!()
    copy_tab!(TAB_XDOT,TAB_YDOT,TAB_XDOT_7,TAB_YDOT_7)
    ####################
    # Step 8
    ####################
    combine_tab!((a70_RK6,a71_RK6,a72_RK6,a73_RK6,a74_RK6,a76_RK6),
                 (TAB_XDOT_1,TAB_XDOT_2,TAB_XDOT_3,TAB_XDOT_4,TAB_XDOT_5,TAB_XDOT_7),
                 (TAB_YDOT_1,TAB_YDOT_2,TAB_YDOT_3,TAB_YDOT_4,TAB_YDOT_5,TAB_YDOT_7),
                  TAB_XDOT,
                  TAB_YDOT)
    drift!(TAB_X_1,TAB_Y_1,
           TAB_X,TAB_Y,
           TAB_XDOT,TAB_YDOT,
           DT)
    TAB_DOT!()
    copy_tab!(TAB_XDOT,TAB_YDOT,TAB_XDOT_8,TAB_YDOT_8)
    ####################
    # Final stage
    ####################
    combine_tab!((b1_RK6,b3_RK6,b4_RK6,b5_RK6,b7_RK6,b8_RK6),
                 (TAB_XDOT_1,TAB_XDOT_3,TAB_XDOT_4,TAB_XDOT_5,TAB_XDOT_7,TAB_XDOT_8),
                 (TAB_YDOT_1,TAB_YDOT_3,TAB_YDOT_4,TAB_YDOT_5,TAB_YDOT_7,TAB_YDOT_8),
                  TAB_XDOT,
                  TAB_YDOT)
    drift!(TAB_X_1,TAB_Y_1,
           TAB_X,TAB_Y,
           TAB_XDOT,TAB_YDOT,
           DT)
    #####
    TIME[1] += 1 # Updating the time
    #####
end