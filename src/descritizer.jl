"""
    descritize(data=[]; bins=[], missing=[])

Descritize a continuous numeric variable to ordinal values {1, 2, ..., N}.

    Inputs:
        data        a vector of continuous values    []
        bins        first value next bin             []
        omit        change these values to missing   []

Example: descritize systolic BP data into three states normal, hypertensive (>= 140) and extremely hypertensive (>= 200) and change entries of 0, -1, and -9 to missing.

descritize(data=myarray, bins=[140, 200], omit=[0, -1, -9])

"""
function descritize(data::Union{Vector,Array}; bins::Vector{<:Number}=[], omit::Array=[])

    low, high = extrema(data)
    high = high + 1
    push!(bins, high)

    pushfirst!(bins, low)
    bc = Int64.(zeros( 1, (length(bins) - 1) ))
    recoded=Array{Union{Int64,Missing}, 1}(undef, size(data,1))
    k = 1
    c = 0
    m =Int64[]
    dl = length(data)
    
    for i in 1:length(bins) - 1
      
        lc = bins[k]
        uc = bins[k+1]

        for d in eachindex(data)
            
            if in(data[d], omit)             
                recoded[d] =  missing
                push!(m, d)
                
            elseif lc <= data[d] < uc

                recoded[d] = k
                c = c + 1
                bc[k] = bc[k] + 1

            end

        end
        
        k = k+1

    end

    mc = length(unique(m))
    tc = mc + c
    println("Input data: $dl\nMissing: $mc\nRecoded: $c\nRecoded counts by bin: $bc\nTotal: $tc")
    
    if tc != dl
        printstyled("WARNING: Not all element properly recoded!\n", color=:yellow)
    end

    return recoded

end

