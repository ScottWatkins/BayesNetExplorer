"""
    plot_network(adjM::Matrix, headerfile::string="", fnode="", gnodes=[], DAG=false, nodeshape=:hexagon )

Plot a Bayesean network from a bnstruct DAG using a DAG adjacency matrix. Node names are taken from the bnstruct header file.

Options:
    gnodes        color nodes given in an array
    fnode         color a feature/target node given as a string 
    method        graph layout [:circular|:stress]
    nodeshape     [:hexagon|:rect|:circle]
"""
function plot_network(adjM::Matrix, headerfile::String=""; fnode="", gnodes=[], DAG=true, nodeshape=:hexagon, method=:stress)

    nodenames = split(replace(readline(headerfile), "\"" =>""), r"\s+")
    
    if length(gnodes) > 0 || length(fnode) > 0
        
        c = []
        ncolors = [:goldenrod1, :cadetblue1, :grey80]
        tn = findall(in([fnode]), nodenames)
        cn = findall(in(gnodes), nodenames)

        if length(fnode) > 0 && length(tn) == 0
            error("\n\nInput fnode ($fnode) not found in network.\n\n")
        end
        
        if length(gnodes) > 0 && length(cn) == 0
            error("\n\nInput gnodes name(s) not found in network.\n\n")
        end

        for i in eachindex(nodenames)
            if in(i, tn)
                push!(c, ncolors[1])
            elseif in(i, cn)
                push!(c, ncolors[2])
            else
                push!(c, :grey80)
            end                        
        end
    else
        c = repeat([:grey80], length(nodenames))
    end

    function cpad(nodenames::Array=[])
        #nodenames = replace.(nodenames, "_" => "")
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
    sz = maximum(length.(nodenames))
    n = length(nodenames)
    w = repeat([sz], n)
    s = repeat([nodeshape], n)
    plot(size(1000,1000))
    
    function make_sym(T)   #convert the bnstruct dag matrix to symetrical
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
              linewidth=1,
              nodecolor=c,
              nodeshape=nodeshape,
              nodesize=0.1,
              #node_weights=w,
              #color=:green,
              curves=false,
              axis_buffer=0.2,
              linecolor = "blue",
              linealpha = 0.99,
              dims=2,
              names=nodenames,
              fontsize=5,
              shorten=0.0,
              method=method
              )
end
