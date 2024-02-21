"""
    benjhoc(pval_vec; fdr=0.05, verbose=true, digits=4)

Calculate the Benjamini-Hockberg corrected p-values for a vector of raw p-values and return the corrected p-values in the input order.
 
P-values < (i/m)*Q are significant, where i is the rank, m is the number of tests, and Q is the false discovery rate.

# Options:
    fdr        False Discovery Rate        0.05
    verbose    Print the ordered table     true
    digits     round to this many digits      4

Reference: McDonald, J.H., Handbook of Biological Statistics

"""
function benjhoc(A=[]; fdr::Float64=0.05, digits::Int64=10, verbose::Bool=true)

    Q = fdr
    A = vec(A)

    An = Symbol.(Tuple([1:1:length(A)...]))
    Ao = NamedTuple{An}(A) #original order  
    
    i = 0
    m = length(A)
    A = sort(A)
    cutoff = Int[]
    bhcorrected = Float64[]
    bhpvalues = Float64[]

    for t in A    
       
        i += 1 
        n = (i/m) * Q
        push!(bhcorrected, n)

        bhp = round(t * m / i, digits=digits)
        push!(bhpvalues, bhp)

        if t < n
            push!(cutoff, i)
        end      

    end

    length(cutoff) == 0 ? push!(cutoff, 10 ,0) : nothing
    
    for e in length(bhpvalues):-1:2
        #adjust bhp taking min values
        ee = min(bhpvalues[e], bhpvalues[e-1])
        bhpvalues[e-1] = ee
    end
    
    line = "-"^50

    if verbose == true
        println(line)
        println("pvalue (i/m)*Q\tsignificance\t B-H p-value")
        println(line)

        for idx in 1:length(bhcorrected)

            bhp = bhpvalues[idx]

            if A[idx] < bhcorrected[idx]
                
                bhc = round(bhcorrected[idx], digits=digits)       
                
                print(A[idx], " < ", bhc,)
       
                if idx ≤ cutoff[end]
                    println("\tsignificant\t$bhp")
                else
                    println("\tnon-significant\t$bhp")
                end
            
            elseif A[idx] ≥ bhcorrected[idx]

                bhc = round(bhcorrected[idx], digits=digits)       
     
                print(A[idx], " ≥ ", bhc)
    
                if idx ≤ cutoff[end]
                    println("\tsignificant\t$bhp")
                else
                    println("\tnon-significant\t$bhp")
                end        

            end
            
        end
        println(line)
    end
    
    c = Dict(zip(A, bhpvalues))
    rv = []
    
    for i in eachindex(Ao) # rv is original order vector 
        v = Ao[Symbol(i)]  # of new pvalues
        np = c[v]
        push!(rv, np)
    end

    return rv
    
end

