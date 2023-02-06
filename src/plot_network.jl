"""
    plot_network(adjM::Matrix, headerfile::String; fnode::String="", gnodes::Array=[], DAG::Bool=true, nodeshape::Symbol=:circle )

    Plot a network using an adjacency matrix. Node names are taken from the bnstruct header file (e.g., BN.header). Required: adjacency matrix and the header file.

    Options
        fnode         color target node           "nodename" 
        gnodes        color nodes in array     ["node1", "node2"]
        bnodes        color nodes in array     ["node1", "node2"]
        cnodes        color node in array      ["node1", "node2"]
        method        graph layout                [:circular|:stress|:hexagon]
        nodeshape     node shape                  [:hexagon|:rect|:circle]
        nodesize      node scalar                 [0.05]
        trimnames     max characters              [20]
        DAG           network is DAG              [true]
        ncolors       node colors; [f,g,b,other]  [:orange, :blue, ...]
"""
function plot_network(adjM::Matrix, headerfile::String=""; fnode="", gnodes=[], bnodes=[], cnodes=[],  DAG=true, nodeshape=:rect, method=:stress, trimnames=10, fontsize::Number=5, nodesize::Float64=0.1, curves::Bool=false, dims::Int64=2, ncolors=[:orange, :skyblue1, :lightgoldenrod1, :palegreen, :grey80], numlabels=false)
    
    nodenames = split(replace(readline(headerfile), "\"" =>""), r"\s+")
    
    if length(gnodes) > 0 || length(fnode) > 0
        
        c = []
        ncolors = ncolors

        tn = findall(in([fnode]), nodenames)
        cn = findall(in(gnodes), nodenames)
        cn2 = findall(in(bnodes), nodenames)
        cn3 = findall(in(cnodes), nodenames)
        
        if length(fnode) > 0 && length(tn) == 0
            error("\n\nInput fnode ($fnode) not found in network.\n\n")
        end
        
        if length(gnodes) > 0 && length(cn) == 0
            error("\n\nInput gnodes name(s) not found in network.\n\n")
        end

        if length(bnodes) > 0 && length(cn) == 0
            error("\n\nInput bnodes name(s) not found in network.\n\n")
        end

        for i in eachindex(nodenames)
            if in(i, tn)
                push!(c, ncolors[1])
            elseif in(i, cn)
                push!(c, ncolors[2])
            elseif in(i, cn2)
                push!(c, ncolors[3])
            elseif in(i, cn3)
                push!(c, ncolors[4])
            else
                #push!(c, :linen)
                push!(c, :snow2)
            end                        
        end
    else
        c = repeat([:grey90], length(nodenames))
    end

    function cpad(nodenames::Array=[])
        #nodenames = replace.(nodenames, "_" => "")

        if trimnames < 20
            for i in eachindex(nodenames)
                if length(nodenames[i]) > trimnames
                    nodenames[i] = nodenames[i][1:trimnames]
                end
            end
        end
        

        sz = maximum(length.(nodenames))   
        pd = Int.(ceil.(( sz .- length.(nodenames)) ./ 2)) .+ 5

        for i in eachindex(nodenames)

            pr = length(nodenames[i]) + pd[i]
            npr = rpad(nodenames[i], pr)
            pl = length(npr) + pd[i]
            npl = lpad(npr, pl)
            nodenames[i] = npl
        end
        return nodenames
    end

    nodenames = cpad(nodenames)

    if numlabels == true
        nodenames = collect(1:1:length(nodenames))
    end

    sz = maximum(length.(nodenames))
    n = length(nodenames)
    w = repeat([sz], n)
    s = repeat([nodeshape], n)
    plot(size(500,500))

    
    function make_sym(T)   #convert the bnstruct dag matrix to symmetrical

        M = zeros(size(T))

        for i=1:size(T,1)
            for j=1:size(T,2)
                if T[i,j] > 0.0
                    M[i,j] = T[i,j]
                    M[j,i] = T[i,j]
                end
                M[i, i] = 0.0 
            end
        end
        return M
    end

    if DAG == true
        M = adjM
    else
        M = make_sym(adjM)
    end

    graphplot(M,
              linewidth=1.5,
              nodecolor=c,
              nodeshape=nodeshape,
              nodesize=nodesize,
              #node_weights=w,
              #color=:green,
              curves=curves,
              axis_buffer=0.15,
              linecolor = :blue,
              linealpha = 0.99,
              dims=dims,
              names=nodenames,
              fontsize=fontsize,
              shorten=0.0,
              method=method
              )
end
