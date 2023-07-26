"""
    moral_net, edges = dag2undirected(DAG)

A function to convert an asymmetric adjacency matrix
representing a directed acyclic graph (DAG) to a
symmetric matrix representing an undirected moralized
network. The directionality of the arrows in the DAG
are removed and a line marrying any unconnected
parents is added. An edgewidth matrix is also returned
"""

function dag2undirected(DAG::Matrix)
    
    if issymmetric(DAG) || sum(diag(x)) > 0
        error("Input matrix should be a DAG (assymmetric) and diagonal values must be zeros.")
    end
    
    Q = deepcopy(DAG) # prevent modification of input matrix!

    newarcs = []

    for c in 1:size(Q,2)    

        idx = findall(Q[:,c] .== 1)  #cols rep. shared offspring              

        if length(idx) > 1
            for i in 1:length(idx) - 1
                for j in 2:length(idx)
                    if i == j        #diag stays 0.0                         
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
            edgM[i[1],i[2]] = 0.2
            edgM[i[2],i[1]] = 0.2
	end
    end

    M = zeros(size(Q))
    
    for i=1:size(Q,1)
        for j=1:size(Q,2)
            if Q[i,j] > 0.0
                M[i,j] = Q[i,j]
                M[j,i] = Q[i,j]
            end
            M[i, i] = 0.0 
        end
    end

    M = Int.(M)
    
    return M, edgM

end
