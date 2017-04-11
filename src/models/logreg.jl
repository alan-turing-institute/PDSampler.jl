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
immutable LogReg
    n::Int           # number of observations
    p::Int           # number of dimensions
    X::Matrix{Float} # design matrix
    y::Vector{Float} # response
    b::Float         # proxy for N*L, L=Lipschitz constant
    # constructor
    LogReg(X,y,b) = new(size(X,1), size(X,2), X, y, b)
end

# robust log(1+exp(x)) (from JuliaStats/StatsFuns.jl/src/basicfuns.jl)
log1pexp(x::Float) = x < 18.0 ? log1p(exp(x)) : x < 33.3 ? x + exp(-x) : x

# logistic and related functions
logistic(x::Float)        =  1.0 / (1.0 + exp(-x))
loglogistic(x::Float)     = -log1pexp(-x)
gradloglogistic(x::Float) =  logistic(-x)

"""
    loglik(lr, w)

Loglikelihood of a Logistic Regression model `lr` at `w`.
"""
function loglik(lr::LogReg, w::Vector{Float})::Float
    sum( loglogistic.(lr.y.*(lr.X*w)) )
end

# helper function, to be used only here; the big term on the right hand side
# is used in all the gradient computations
glli(w, lr, i) = gradloglogistic(lr.y[i]*(dot(lr.X[i,:],w)))

"""
    gradloglik(lr, w)

Gradient of the loglikelihood of a Logistic Regression model `lr` at `w`
"""
function gradloglik(lr::LogReg, w::Vector{Float})::Vector{Float}
    sum( lr.y[i]* glli(w,lr,i) * lr.X[i,:]  for i in 1:lr.n)
end

"""
    gradloglik_cv(lr, wstar)

Return a function that can computes a gradient of the loglikelihood for the
model `lr` around a point `wstar`.
"""
function gradloglik_cv(lr::LogReg, wstar::Vector{Float})::Function
    # Unbiased estimate of \nabla U_j(x) using control variates
    gll_star = gradloglik(lr, wstar)
    function gll_cv(w::Vector{Float}, i=rand(1:lr.n))::Vector{Float}
        tw, tws = glli(w,lr,i), glli(wstar,lr,i)
        # unbiased estimate of the gradient
        gll_star + lr.n*lr.y[i] * (tw*lr.X[i,:]-tws*lr.X[i,:])
    end
end
