"""
    Dependent function. This function is used within bne.jl.
"""
function ConProb(;f="", g=g, gs=gs, type="conditional", vars=vars, verbose=false, rr_bootstrap=rr_bootstrap) # get P(feature|g1, g2, ... gn)

    if length(g) != length(gs)
        println("g was: ", length(g), " gs was: ", length(gs))
        error("Please define a state for each conditional feature!")
    end

    if !issubset(g, vars)
        error("At least one input values ($g) is not a known variable!")
    end

    q = "\'"

    feature = join(["c(", q, f, q, ")"], "") #gRain string   

    gnodes=""; gstates="";

    for i in g
        gnodes = join([gnodes join([q,i,q], "")], ",")
    end

    gnodes= join(["c(", gnodes[2:end], ")"], "")

    for j in gs
        gstates = join([gstates join([q,j,q], "")], ",")
    end
    
    gstates = join(["c(", gstates[2:end], ")"] , "") # done creating strings 
    R"f <- $feature"
    R"g <- $gnodes"; R"gs <- $gstates"
    R"ev_net1 <- setEvidence(pnet, nodes = $g, states = $gs )"   #set states
    grq = rcopy(R"grq = querygrain(ev_net1, nodes = $f, type=$type)"); #P(feature|g1,g2,...)

    if verbose == true
        println("$('-'^75)")
        println("P($f)|$g\nConditional state(s): $gs\n")
        R"print(grq)"
        println("$('-'^75)")
    end
    
    probout = round.(grq, digits=6)

    qout = join([f, probout, join(string.(g), ","), join(string.(gs), ",") ], "|")

    if rr_bootstrap > 0
        qout = probout
    end

    return qout

end
