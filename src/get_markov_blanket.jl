"""
    get_markov_blanket(DAG, features::Array, target::Union{String,Int})

Return the markov blanket of a target node for a DAG. Features is an array with the column names of the variables in the DAG. The target is the integer or name that specifies the target's col id in the DAG. The markov blanket contains the target's parents and offspring and the parents of the offspring. Input DAG format, D[i,j] must be D[parent,offspring], that is, arrow goes from row to column. Multiple targets, input in an array, are allowed with the output being the union of all blankets. 
"""
function get_markov_blanket(DAG::Matrix, features::Array{String}, target::Union{String,Int, Array} )
    
    function process_target(DAG::Matrix, features::Array{String}, target::Union{String,Int,Array} )

        if  typeof(target) == String
            target = findall(occursin.(Regex("^$target\$"), features))
            target = target[1]
            println("Target index is $target")
        elseif typeof(target) == Vector{Int}
            label = features[target]
            println("Target indices are $target")
        end
        
        mb = Int.(zeros(size(DAG)))  
        
        mb[:, target] = DAG[:, target]  # all parents of target
        mb[target, :] = DAG[target, :]  # all offspring of target
    
        for i in 1:size(DAG, 1)
            for j in 1:size(DAG, 2)
                if i == target && DAG[i,j] > 0 # is offspring
                    if sum(DAG[:,j]) > 1   # offspring has parent(s)
                    mb[:, j] = DAG[:, j]   # add parent of offspring of target 
                    end
                end
            end
        end
        
        mbn = []
    
        for i in 1:size(mb,1)
            if length(findall(mb[i,:] .> 0)) > 0
                push!(mbn,i)
            end
            if length(findall(mb[:,i] .> 0)) > 0
                push!(mbn,i)
        end
        end
        
        mbnames = features[mbn]
        mbnames = sort(unique(mbnames))
        
        return (mb, mbnames)
        
    end

    fmb = Int.(zeros(size(DAG)))
    fmbnames = []

    for i in target     # Iterate over multiple targets
        mb, mbnames = process_target(DAG, features, i )
        fmb = fmb .+ mb
        fmbnames = push!(fmbnames, mbnames)
    end

    return (fmb, fmbnames)

end

