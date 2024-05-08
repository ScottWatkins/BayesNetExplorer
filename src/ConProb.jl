"""
    Dependent function. This function is used within bne.jl to calculate conditional probability. Note user defined rrr denominator is passed in as gs. This function cycles twice in bne to get the initial query and then the negated query for the relative risk denominator.

"""
function ConProb(;f=f, fs=fs, g=g, gs=gs, type=type, vars=vars, verbose=verbose, rr_bootstrap=rr_bootstrap, query=query) # get P(feature|g1, g2, ... gn)

    if length(g) != length(gs)
        println("\nMismatched conditional query:\ng elements: ", length(g), " but gs elements: ", length(gs))
        error("\nPlease match each conditional feature with a conditional state!\n\n")
    end

    if !issubset(g, vars)
        error("\n\nAt least one input values ($g) is not a known variable!\n\n")
    end

    q = "\'"

    fnodes=""; gnodes=""; gstates="";

    if typeof(f) == Vector{String}
        for i in f
            fnodes = join([fnodes join([q,i,q], "")], ",")
         end
        fnodes= join(["c(", fnodes[2:end], ")"], "")
    else 
        fnodes = join(["c(", q, f, q, ")"], "") #gRain string   
    end
    
    for i in g
        gnodes = join([gnodes join([q,i,q], "")], ",")
    end

    gnodes= join(["c(", gnodes[2:end], ")"], "")


    for j in gs
        gstates = join([gstates join([q,j,q], "")], ",")
    end

    gstates = join(["c(", gstates[2:end], ")"] , "") # done creating strings 

    R"f <- $fnodes"
    R"g <- $gnodes"; R"gs <- $gstates"

    # Here the compiled net get set for gstates of gnodes, that is smoke=yes means
    # the net is set for P(smoke=1.0). The net is then propogated.
    # println("ev_net1 <- setEvidence(pnet, nodes = $gnodes, states = $gstates, propagate=TRUE )" )

    R"ev_net1 <- setEvidence(pnet, nodes = $g, states = $gs, propagate=TRUE )"   #set states

    # The final query is made for one or more target nodes
    R"querygrain(ev_net1, nodes = $f, type=$type)"

    R"grq = querygrain(ev_net1, nodes = $f, type=$type)"
    dname = rcopy(R"dname = dimnames(grq)")  #dname is ordered dict from grain

    ff = Symbol.(f)

    if typeof(ff) == Symbol 
        ffs = Dict(ff => fs)        
    elseif typeof(ff) == Vector{Symbol}
        ffs = Dict(zip(ff, fs))    # input targets
    end
    
    gr_idx = []                # empty index
                               # Input target-state must be *matched* to grain matrix!
                               # population new index based on grain

    if type == "joint" || type == "conditional"
        for k in keys(dname) 
            i = get(ffs, k, 0)
            if typeof(i) == String
                i = parse(Int64, i)
            end
            push!(gr_idx, i)

        end
    end
    
    grq = rcopy(R"grq = querygrain(ev_net1, nodes = $f, type=$type)");
        
    #Note: ev_net1 is propagated, qrq query P(feature|g1,g2,... is from the propagated net)

    c = ""
    
   if verbose == true && rr_bootstrap == 0 
       println("$('-'^75)")
       for i in 1:length(g)
           c = c * g[i] * "=" * gs[i] * ","
       end
       c = c[1:end-1]
       
       println("$query\n")
       println("Listing all target probabilities...")
       R"print(grq)"
       println("$('-'^75)")

   end
    
    probout = grq

    jpout = join([f, probout, join(string.(g), ","), join(string.(gs), ",") ], "|")
        
    if rr_bootstrap > 0
        boot_pF = probout[CartesianIndex(Tuple(gr_idx))]   #index to grain target probability
        probout = boot_pF
        jpout = probout
    end
    
    return probout, jpout, gr_idx

end
