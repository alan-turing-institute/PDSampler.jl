using PDMP, Base.Test

# ------------------------------------------------------------------------------
# Quick version of functions (for testing) (drawn from some of SJV's old code)
q_loglik(lr,w)     = sum(-log(1+exp(-lr.y.*(lr.X*w))))
q_gradloglik(lr,w) = sum( lr.y[i]*(1/(1+exp(lr.y[i]*dot(lr.X[i,:],w))))*X[i,:]
                            for i in 1:lr.n )
function q_gradloglik_cv(lr,wstar)
    gll_star = q_gradloglik(lr,wstar)
    function gll_cv(w::Vector{Float};i=rand(1:lr.n))
        gll_star+lr.n*(
            lr.y[i] .* 1/(1+exp(lr.y[i] .* (dot(lr.X[i,:],w)))) .* lr.X[i,:] - lr.y[i] .* 1/(1+exp(lr.y[i] .* (dot(lr.X[i,:],wstar)))).*lr.X[i,:])
    end
end
# ------------------------------------------------------------------------------
srand(123)

xvals = 10000*randn(500)

safemaximum(x,s) = length(s)==1?max(x,s):max(x,s[1])

@test maximum(logistic(x)-1.0/(1.0+exp(-x)) for x in xvals) <= 1e-10
vals = [loglogistic(x)-(-log(1.0+exp(-x))) for x in xvals]
@test maximum(vals[!isinf(vals)]) <= 1e-10
@test maximum(norm(gradloglogistic(x)-1.0/(1.0+exp(x))) for x in xvals) <= 1e-10

# Generate data for Likelihood testing etc.
n = 10000           # number of observations
p = 100             # number of dimensions (covariates)
X = randn(n,p)+0.1  # feature matrix
w = 10*rand(p)      # true vector of parameters
# observations according to a logistic thresholded to {-1,1}
y = (logistic.(X*w) .> rand(n)) .* 2.0 .- 1.0
# proxy for N*L upper bound (see PDMP paper TODO document better, check SQUARE?)
b  = sum( mapslices(_->norm(_)^2,X,1) )/4
lr = LogReg(X,y,b)

xvals = [randn(p) for i in 1:5]

@test maximum(loglik(lr,x)-q_loglik(lr,x) for x in xvals) <= 1e-10

@test maximum(norm(gradloglik(lr,x)-q_gradloglik(lr,x)) for x in xvals) <= 1e-10

xstar = randn(p)
i = rand(1:n)

# some aggregation of numerical noise
@test maximum(norm(gradloglik_cv(lr,xstar)(x,i) -
                    q_gradloglik_cv(lr,xstar)(x,i=i)) for x in xvals) <= 1e-10
