export
    LogReg,
    logistic,
    loglogistic,
    gradloglogistic

"""
    LogReg(X, y, b)

Data Model for a Logistic Regression with `X` the design matrix, `y` the
response in `(-1,1)` and `b` is the proxy for the Lipschitz constant in order
to use thinning in the BPS algorithm.
"""
struct LogReg
    n::Int # number of observations
    p::Int # number of dimensions
    X::Matrix{Real} # design matrix
    y::Vector{Real} # response
    b::Real  # proxy for N*L, L=Lipschitz constant
    # constructor
    LogReg(X,y,b) = new(size(X,1), size(X,2), X, y, b)
end

# robust log(1+exp(x)) (from JuliaStats/StatsFuns.jl/src/basicfuns.jl)
log1pexp(x::Real) = x < 18.0 ? log1p(exp(x)) : x < 33.3 ? x + exp(-x) : x

# logistic and related functions
logistic(x::Real) = 1.0 / (1.0 + exp(-x))
loglogistic(x::Real) = -log1pexp(-x)
gradloglogistic(x::Real) = logistic(-x)

"""
    loglik(lr, w)

Loglikelihood of a Logistic Regression model `lr` at `w`.
"""
loglik(lr::LogReg, w::Vector{<:Real}) = sum(loglogistic.(lr.y .* (lr.X * w)))

# helper function, to be used only here; the big term on the right hand side
# is used in all the gradient computations
glli(w, lr, i) = gradloglogistic(lr.y[i] * (dot(lr.X[i, :], w)))

"""
    gradloglik(lr, w)

Gradient of the loglikelihood of a Logistic Regression model `lr` at `w`
"""
gradloglik(lr::LogReg, w::Vector{<:Real}) =
    sum(lr.y[i] * glli(w,lr,i) * lr.X[i,:] for i ∈ 1:lr.n)

"""
    gradloglik_cv(lr, wstar)

Return a function that can computes a gradient of the loglikelihood for the
model `lr` around a point `wstar`.
"""
function gradloglik_cv(lr::LogReg, wstar::Vector{<:Real})
    # Unbiased estimate of ∇U_j(x) using control variates
    gll_star = gradloglik(lr, wstar)

    function gll_cv(w::Vector{<:Real}, i=rand(1:lr.n))
        tw, tws = glli(w, lr, i), glli(wstar, lr, i)
        # unbiased estimate of the gradient
        gll_star + lr.n * lr.y[i] * (tw * lr.X[i, :] - tws * lr.X[i,:])
    end
    return gll_cv
end
