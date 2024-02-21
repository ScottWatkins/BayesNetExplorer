"""
    tdistCI(A, alpha=0.95)

Get the upper and lower CI values for an input array using a T-distribution
"""
function tdistCI95(A; alpha::Float64=0.95)

    x = mean(A)        #sample mean
    σ = std(A)         #sample standard deviation
    n = size(A, 1)     #sample size
    ν = n - 1          #degrees of freedom
    α = alpha          #
    
    ζ =  quantile(TDist(ν), 1.0 - α/2.0) #critical value at alpha
    println(ζ)
    (CI_L, CI_U) = (x - ζ * σ/√n, x + ζ * σ/√n)  #lower and upper value in A

    return CI_L, CI_U
    
end
