"""
    Dependent function. This function is used within bne.jl to calculate conditional probability. Note user defined rrr denominator is passed in as gs. This function cycles twice in bne to get the initial query and then the negated query for the relative risk denominator.

"""
function ConProb(;f=f, fs=fs, g=g, gs=gs, type=type, vars=vars, verbose=false, rr_bootstrap=rr_bootstrap, query=query) # get P(feature|g1, g2, ... gn)

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

    # NOTES: setting type="marginal" for P(xray=2,bronc=2|smoke=2) returns 
    # values similar values as cpq for P(xray=1,2|smoke=2) and P(bronc=1,2|smoke=2).
    # For "P(xray=2|smoke=2,bronc=2)", marginal and conditional match cpq 
    # For "P(A=2,B=1,C=2,D=1|Y=2,Z=2)" the final P-value is in grq[2,1,2,1], etc.
    # The propagated values can be very different than cpq especially if
    # ANY nodes are not in same markov blanket.  
    
    #println("f $f | fnodes $fnodes | g $g |gnodes: $gnodes | gs $gs |  gstates: $gstates")
    
    R"f <- $fnodes"
    R"g <- $gnodes"; R"gs <- $gstates"

    # Here the compiled net get set for gstates of gnodes, that is smoke=yes means
    # the net is set for smoke=1.0. The net is then propogated.
    #println("ev_net1 <- setEvidence(pnet, nodes = $gnodes, states = $gstates, propagate=TRUE )" )

    R"ev_net1 <- setEvidence(pnet, nodes = $g, states = $gs, propagate=TRUE )"   #set states

    # The final query is made for one or more target nodes
    grq = rcopy(R"grq = querygrain(ev_net1, nodes = $f, type=$type)");

    #Note: ev_net1 is propagated, qrq query P(feature|g1,g2,... is from the propagated net)

    c = ""
    
    if verbose == true
        println("$('-'^75)")
        for i in 1:length(g)
            c = c * g[i] * "=" * gs[i] * ","
        end
        c = c[1:end-1]

        
        println("$query\n")
        R"print(grq)"
        println("$('-'^75)")
    end
    
    probout = grq
    
    jpout = join([f, probout, join(string.(g), ","), join(string.(gs), ",") ], "|")

    if rr_bootstrap > 0
        boot_pF = probout[CartesianIndex(Tuple(fs))] #index to allow multi-target queries
        probout = boot_pF
        jpout = probout
    end
    
    return probout, jpout

end
