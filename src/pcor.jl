"""
    pmat, pairs = pcor(cm::Union{Matrix,Float64}, obs; pval=0.05, features=[])

Calculate the two-tailed p-values for a correlation matrix based on n observations. Get (r, p-value) pairs significant at a given p-value threshold. P-values are calculated on a T-distribution. Input may also be a single correlation.
"""
function pcor(cm::Union{Matrix,Float64}, obs::Int64; pval::Float64=0.05, features::Array=[])

    if typeof(cm) == Float64   #process single input correlation
        cm = [1.0 cm; cm 1.0]
    end
    
    
    if size(cm,1) !== size(cm,2)
        error("Input matrix is not square. Input must be a square correlation matrix.")
    end

    if obs == 0
        error("You must enter the number of observations (obs) used to calculate the correlation")
    end
    
    obs == 2 ? n = obs + 1 : n = obs

    pmat= zeros(size(cm))
    td = TDist(n-2)
    acm = abs.(cm)

    for i in 1:size(acm, 1)
        for j in 1:size(acm, 2) 
            t = acm[i,j] * sqrt( (n - 2) / (1 - acm[i,j]^2) )
            p = ccdf(td,t) * 2
            pmat[i,j] = p
        end
    end
    
    sigidx = findall(tril(pmat .<= pval, -1))

    if length(features) == size(cm,1)
        idxvals1 = [ features[ i[1] ] for i in sigidx]
        idxvals2 = [ features[ i[2] ] for i in sigidx]
    else
        idxvals1 = [ i[1] for i in sigidx]
        idxvals2 = [ i[2] for i in sigidx]
    end
    
    cmvals = [ cm[i[1], i[2]] for i in sigidx]
    pvals = [ pmat[i[1], i[2]] for i in sigidx]
    sighits = hcat(idxvals1, idxvals2, cmvals, pvals)

    return pmat, sighits

end



