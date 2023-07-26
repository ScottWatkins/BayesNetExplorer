"""
    Creates a line marrying unconnected parents.
    The input matrix should be a DAG.


"""
function moralizeDAG(adjM; moral_line::Float64=1.0)

    Q = deepcopy(adjM)

    newarcs = []

    for c in 1:size(Q,2)

        idx = findall(Q[:,c] .== 1)  #cols rep. shared offspring              

        if length(idx) > 1
            for i in 1:length(idx) - 1
                for j in 2:length(idx)
                    if i == j        #diag stays zero      
                        continue
                    end
                    push!(newarcs, [idx[i], idx[j]])
                end
            end
        else
            continue
        end
    end

    edgM = rand([1.5], size(Q))
    
    for i in newarcs
        if Q[i[1],i[2]] == 0 && Q[i[2],i[1]] == 0    
            Q[i[1],i[2]] = 1
            Q[i[2],i[1]] = 1
            edgM[i[1],i[2]] = moral_line
            edgM[i[2],i[1]] = moral_line
        end
    end

    return Q, edgM

end
