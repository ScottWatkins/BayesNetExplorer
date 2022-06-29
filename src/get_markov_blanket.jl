"""
    get_markov_blanket(DAG, features::Array, target::Union{String,Int})

Return the markov blanket of a target node for a DAG. Features is an array with the column names of the variables in the DAG. The target is the integer or name that specifies the target's col id in the DAG. The markov blanket contains the target's parents and offspring and the parents of the offspring. Input DAG format, D[i,j] must be D[parent,offspring].
"""
function get_markov_blanket(DAG::Matrix, features::Array{String}, target::Union{String,Int}, )
    
    if typeof(target) == String
        target = findall(occursin.(target, features))
        target = target[1] 
        println("Target index is $target")
    elseif typeof(target) == Int
        label = features[target]
        println("Target label is $label")
    end
    mb = Int.(zeros(size(DAG)))  

    mb[:, target] = DAG[:, target]  # all parents of target
    mb[target, :] = DAG[target, :]  # all offspring of target
    
    for i in 1:size(DAG, 1)
        for j in 1:size(DAG, 2)
            if i == target && DAG[i,j] > 0 # is offspring
                if sum(DAG[:,j]) > 1       # offspring has parent(s)
                    mb[:, j] = DAG[:, j]   # add parent of offspring of target 
                    #println("==>", i, " ", j, " ", DAG[:,j])
                end
            end
        end
    end
    
    return (mb)
    
end

