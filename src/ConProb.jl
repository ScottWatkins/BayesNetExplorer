"""
    Dependent function. This function is used within bne.jl to calculate conditional probability. Note user defined rrr denominator is passed in as gs.
"""
function ConProb(;f="", fs=fs, g=g, gs=gs, type="conditional", vars=vars, verbose=false, rr_bootstrap=rr_bootstrap) # get P(feature|g1, g2, ... gn)

    if length(g) != length(gs)
        println("\nMismatched conditional query:\ng elements: ", length(g), " but gs elements: ", length(gs))
        error("\nPlease match each conditional feature with a conditional state!\n\n")
    end

    if !issubset(g, vars)
        error("\n\nAt least one input values ($g) is not a known variable!\n\n")
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
    R"ev_net1 <- setEvidence(pnet, nodes = $g, states = $gs, propagate=TRUE )"   #set states
    grq = rcopy(R"grq = querygrain(ev_net1, nodes = $f, type=$type)");

    #Note: ev_net1 is propagated and qrq query P(feature|g1,g2,... is from the propagated net)
    c = ""
    if verbose == true
        println("$('-'^75)")
        for i in 1:length(g)
            c = c * g[i] * "=" * gs[i] * ","
        end
        c = c[1:end-1]
        
        println("P($f=$fs|$c)\n")
        R"print(grq)"
        println("$('-'^75)")
    end
    
    probout = round.(grq, digits=8)

    jpout = join([f, probout, join(string.(g), ","), join(string.(gs), ",") ], "|")

    if rr_bootstrap > 0
        jpout = probout
    end

    return probout, jpout

end
