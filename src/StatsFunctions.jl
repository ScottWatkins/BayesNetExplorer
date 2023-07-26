#module StatsFunctions

# Functions included in this module
# 1. Normalization functions - box-cox, min-max, decimal
# 2. P-value for Pearson correlation
# 3. Benjamini-Hochberg correction for multiple comparisions
# 4. Generalized Z test for proportions in a single sample

using Distributions

#export norm_dat, cor_p, bh_correction, OneSampleZtest

#----------------------------------------------------------
# Normalization functions
#----------------------------------------------------------
"""
    norm_dat(x::Array{Number}, method=:minmax)

Normalize an array of numbers X so that all elements of X
are between 0 and 1, or -1 and 1, etc. using standard
normalization functions.

Available methods

    :minmax  [default]

        Normalized range (0 to 1)
        X[1..n] = (X[n] - min[X])/(max[X] - min[X])

    :boxcox, lamda=0.5  (square root)

        Power transformation
        X[1..n] = (X[n]^lambda - 1) / lamda

        Specify lamda (-2 through 2)
        e.g.  0.5 => square root
             -0.5 => inverse square root 
                2 => square
               -2 => inverse square
                   
    :decimal

        Normalized range (-1 to 1)
        X[1..n] = X[n] / 10^(max. num. X[n] digits)


The order of the array is not altered.

"""
function norm_dat(x::Array{<:Number}; method=:minmax, lamda=0.5)

    min,max = extrema(x)

    if  method === :boxcox

        norm = ((x.^lamda) .- 1) ./ lamda

        return norm

    elseif method === :decimal

        if length(digits(max)) > length(digits(min))
            d = length(digits(max))
        else
            d = length(digits(min))
        end

        norm = (x ./ 10^d)

        return norm

    elseif method === :minmax

        d = max-min
        norm = (x .- min) ./ d
        return norm

    end    

end
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------



#----------------------------------------------------------
# P-values for Pearson correlations
#----------------------------------------------------------

"""                                       
    cor_p(r, n, t)
    
Get the p-value for a Pearson correlation.
P is based on a T-distribution.
r is the correlation, n is number of values,
and t is number of tails (2,1).
"""                                                                                 
function cor_p(r, n, t)

    r = abs.(r)
    df = n - 2

    z = (r*sqrt(n-2))/(sqrt(1-r^2))
    p = round.(ccdf(TDist(df), z), digits=4)

    if t == 2
        p = p * 2
    end

    return(p)
    
    # emperical, slow
    #  x=rand(TDist(df), 100000)
    #  round.(((length(findall(x .> z)) * t)) /100000, digits=4)
  
end

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------


#----------------------------------------------------------
# Benjamini-Hochberg correction for multiple comparisions
#----------------------------------------------------------

"""
    bh_correction(A, fdr)

 Calculate the Benjamini-Hockberg corrected p-values and show
 the significance for a set of observed p-values given 
 in array A. User sets the false discovery rate (fdr) [0.05].
 
 P-values < (i/m)*Q are significant, where i is the rank, 
 where m is the number of tests, and Q is the false discovery rate.

 All values < the maximum significant p-value are considered
 significant for a given fdr.  

 Reference: McDonald, J.H., Handbook of Biological Statistics
"""
function bh_correction(A, fdr)

    @isdefined(fdr) ? Q = fdr : fdr = 0.05

    A = vec(A)
       
    i = 0
    m = length(A)
    sort!(A)
    cutoff = Int[]
    bhcorrected = Float64[]
    bhpvalues = Float64[]

    for t in A    
       
        i += 1 
        n = (i/m) * Q
        push!(bhcorrected, n)

        bhp = round(t * m / i, digits = 8)
        push!(bhpvalues, bhp)

            if t < n
                push!(cutoff, i)
            end      

    end

    for e in length(bhpvalues):-1:2           #adjust bhp taking min values
        ee = min(bhpvalues[e], bhpvalues[e-1])
        bhpvalues[e-1] = ee
    end

    a = "-"^50
    println(a)
    println("pvalue (i/m)*Q\tsignificance\t B-H pvalue")
    println(a)

    for idx in 1:length(bhcorrected)

        bhp = bhpvalues[idx]

        if A[idx] < bhcorrected[idx]
            
            bhc = round(bhcorrected[idx], digits=4)       

            print(A[idx], " < ", bhc,)
       
            if idx ≤ cutoff[end]
                println("\tsignificant\t$bhp")
            else
                println("\tnon-significant\t$bhp")
            end
               
        elseif A[idx] ≥ bhcorrected[idx]

            bhc = round(bhcorrected[idx], digits=4)       
     
            print(A[idx], " ≥ ", bhc)
    
            if idx ≤ cutoff[end]
                println("\tsignificant\t$bhp")
            else
                println("\tnon-significant\t$bhp")
            end        

        end

    end

end
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------


#----------------------------------------------------------
#--------Generalized Z-score test of proportions-----------
#----------------------------------------------------------
"""
    OneSampleZtest(obs, exp, num)

Calculate the z-score and 2-tailed p-value for an observed
percent compared to its expected percent, where obs is the 
observed fraction, exp is the expected fraction, and num
is the number of samples. The data must be normally
distributed.

Example:

(z, p) = OneSampleZtest(0.59388, 0.50000, 109)

(1.9603, 0.05)

"""
function OneSampleZtest(obs::Float64, exp::Float64, num::Int)

    p::Float64  = obs
    pexp::Float64 = exp
    n::Int64 = num
    alpha = 1.96
    
    se = sqrt( ((pexp * (1 - pexp))) / n )

    z = round((p - pexp) / se, digits=4)

    if z < 0
        pval = round(2 * ( (cdf(Normal(), z))), digits=4)
    else
        pval = round(2 * (1 - (cdf(Normal(), z))), digits=4)
    end

    CIhigh  = round(p + (alpha * sqrt( (p * (1-p)) / n  )), digits=4)
    CIlow = round(p - (alpha * sqrt( (p * (1-p)) / n  )), digits=4)


    println("Z-score = $z")
    println("P-value = $pval")
    println("CI95%: $CIlow - $CIhigh")
    if CIlow <= pexp <= CIhigh
        println("Fail to reject null hypothesis: the 95% CI overlaps $pexp\n")
    else
        println("Reject null hypothesis: the 95% CI does not overlap $pexp.\n")
    end
    
    return (z, pval)

end


"""
    TwoSampleZtest(obs1, obs2, num1, num2)
Calculate the z-score and 2-tailed p-value for the observed
count of sample1 to the observed counts of sample2 given 
num1 number of counts of sample1 and num2 counts of sample2.

Example:
(z, p) = TwoSampleZtest()  
(1.9603, 0.05)
"""
function TwoSampleZtest(obs1, obs2, num1, num2)
    println("Code in progress")
end

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------






#end  #end module
