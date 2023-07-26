"""
    plot_network(adjM::Matrix; headerfile::String="", fnode::String="", gnodes::Array=[], DAG::Bool=true, nodeshape::Symbol=:circle )

    Plot a network using an adjacency matrix. Node names can be taken from the first line of a space delimited file (e.g., BN.header) as the. Required: adjacency matrix.

    Options
        headerfile    space delim file            "filename"
        fnode         color target node           "nodename" 
        gnodes        color nodes in array        ["node1", "node2"]
        bnodes        color nodes in array        ["node1", "node2"]
        cnodes        color node in array         ["node1", "node2"]
        method        graph layout                [:stress|:sfdp|:circular|
                                                  :shell|:spectral|:spring|
                                                  :tree|:buchheim|:arcdiagram|
                                                  :chorddiagram]
        nodeshape     node shape                  [:hexagon|:rect|:circle|
                                                  :ellipse]
        trimnames     max characters              20
        DAG           network is DAG              true
        moralize      join parents if DAG=false   false
        ncolors       node colors [f,g,b,other]   [:orange, :blue, ...]
        nodeweight    vec of relative node sizes  [1.0 ...]
        nodesize      node scalar                 [0.05]
        nodenames     optional node names         []
                      (input must be string)
        linewidth     size of lines               1.5
        curves        curve the arcs              false
        curvescale    amplitude of curve          0.03
        dims          2 or 3-dimensions           2
        edgelabel     matrix for edge labels      []
        edgewidth     matrix of edge widths       [2.0 ...; ... 2.0]               
        boxcolor      color for edge label box    :white
        moral_line    width of add moral lines    1.0

Notes:
1. Edgelabel is a string matrix that matches the dag adjM.
   1. "" indicates no edge label.
   2. If the adjM position has a 1, the edgelabel matrix
      can have a corresponding label. See (2) for input info.
   3. The edgelabels matrix will be automatically converted
      to symmetric if DAG=false.
2. Input non-symetric matrices have directional arcs.
      From row index -> to col index
3. If the matrix is symmetrical, the graph is non-directional.      

"""
function plot_network(adjM::Matrix; headerfile::String="", fnode::String="", gnodes::Array=[], bnodes::Array=[], cnodes::Array=[],  DAG::Bool=true,  nodeshape::Symbol=:circle, method::Symbol=:stress, trimnames::Int=10, fontsize::Number=5, nodesize::Float64=0.06, curves::Bool=false, dims::Int64=2, ncolors::Array=[:orange, :skyblue1, :lightgoldenrod1, :palegreen3, :grey80], moralize::Bool=false, nodenames::Array=[], curvescale::Float64=0.03, edgelabel::Array=[], edgewidth::Array=[], node_weights::Vector=[], boxcolor::Symbol=:white, linewidth::Float64=1.0, moral_line::Float64=1.0 )
    

    if isfile(headerfile)
        nodenames = split(replace(readline(headerfile), "\"" =>""), r"\s+")
    elseif length(nodenames) == size(adjM,2)
        println("Using user specified nodenames...")
        println(nodenames)
    else
        println("Using sequential numbering for nodes names...")
        nodenames = string.(collect(1:1:size(adjM,2)))
        println(nodenames)
    end

    if DAG == true && moralize == true
        error("\nMoralization converts a directed network to a non-directed network,\nso DAG should be set to false.\n")
    end

    
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

        if length(bnodes) > 0 && length(cn2) == 0
            error("\n\nInput bnodes name(s) not found in network.\n\n")
        end

        if length(cnodes) > 0 && length(cn3) == 0
            error("\n\nInput cnodes name(s) not found in network.\n\n")
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

        M = Int.(M)

        return M
    end
    
    if length(edgewidth) == 0         # initialize edgewidth matrix
        edgewidth = rand([linewidth], size(adjM)) 
    end

    if length(node_weights) == 0      # initialize node_weight vector     
        node_weights = rand([1.0], size(adjM,1)) 
    end

    edgM = rand([1.5], size(adjM))
    
    if DAG == true

        M = adjM

    else

        if moralize == true
            
            adjM_moral, edgM = moralizeDAG(adjM, moral_line=moral_line)  #marry parents        
            printstyled("Creating non-directional, moralized network...\n", color=:green)
            M = make_sym(adjM_moral)       #func to make symmetrical

        else

            M = make_sym(adjM)
            printstyled("INFO: Creating non-directional but non-moralized network...\n", color=:green)

        end

    end
    
    graphplot(M,
              linewidth=linewidth,
              names=nodenames,
              nodecolor=c,
              nodeshape=nodeshape,
              nodesize=nodesize,
              curves=curves,
              curvature_scalar=curvescale,
              axis_buffer=0.1,
              linealpha = 0.99,
              dims=dims,
              fontsize=fontsize,
              shorten=0.0,
              edgelabel=edgelabel,
              edgelabel_offset=0.0,
              edge_label_box=true,
              self_edge_size=0.1,
              linecolor=boxcolor,    #for edge label box 
              edgewidth=edgM,
              node_weights=node_weights,
              arrow=:closed,
              method=method
              )

end
