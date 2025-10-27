##################################################
# Function used to find zero
# It is a simple bisection algorithm, but it makes no allocations and is sufficiently fast
# It allows us not to have to use the Roots library that makes a lot of allocations
# The optional tolerances are set to the same as the ones I found in Roots.jl.
# @ATTENTION, the tolerances are most likely overkill -- it may prevent convergence for high-order resonances
# @IMPROVE -- it could be a good idea to put a counter of iterations, to prevent the algorithm from getting stuck?
##################################################
function bisection(fun, xl::Float64, xu::Float64, tolx::Float64=100.0*eps(Float64), tolf::Float64=100.0*eps(Float64))
    if (xl > xu)
        xl, xu = xu, xl # Ordering the input arguments
    end
    #####
    fl, fu = fun(xl), fun(xu) # Evaluating the function on the bracket
    #####
    if (abs(fl) <= tolf) # We have already found a solution on the left bracket
        return xl # Returning the left bracket
    end
    #####
    if (abs(fu) <= tolf) # We have already found a solution on the right bracket
        return xu # Returning the right bracket
    end
    #####
    @assert fl*fu < 0.0 "bisection: NOT A BRACKET"
    #####
    while true # Bisection loop
        #####
        xm = (xl+xu)*0.5 # Middle value
        #####
        if (abs(xu-xl) <= tolx) # The considered bracket is smaller than the tolerance
            return xm # Returning the middle value
        end
        #####
        fm = fun(xm) # Value at the midpoint
        #####
        if (abs(fm) <= tolf) # The middle value is below the threshold
            return xm # Returning the middle value
        end
        #####
        # Otherwise, we iterate the bisection
        if (fm*fl < 0.0) # One root can be found between the left point and the middle point
            xu, fu = xm, fm # The upper point becomes the midpoint
        else
            xl, fl = xm, fm # The lower point becomes the midpoint
        end
    end
end