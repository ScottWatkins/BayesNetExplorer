"""

    bootstrapRRtable(RRtable)

Bootstrap the output of the RRcalculator using bne.

"""
function bootstrapRRtable(RRtable::Union{String,DataFrame}; rr_bootstrap::Int=100, bootout::String="bootstrap_results.csv")

    #addprocs(threads)

    if isfile(bootout)
        error("\n\nFound existing file call $bootout.\nPlease remove or use a different filename.\n\n")
    else
        OUT = open(bootout, "w")
        println(OUT, "Target\tConditionals\tP(T)\tP(T|C)\tARR\tARRest\tARRlower\tARRupper\tRRR\tRRRest\tRRRlower\tRRRupper")
        close(OUT)
    end
                
    line="-"^70
    
    if typeof(RRtable) == DataFrame
        df = RRtable
    elseif isfile(RRtable)
        df = CSV.read(RRtable, DataFrame, types=String)
    end

    q_all = querywriter(df)

    f_all = df[!,1]

    fs = string(split(names(df)[2], r"_|\|" )[2]  )

    if fs == "Yes"
        println("Setting fs=Yes to fs=2...")
        fs = "2"
    end
    
    g_all = df[!,9]
    gs_all = df[!,11]

    for i in eachindex(f_all)

        f = f_all[i]
        fs = fs
        g = string.(split(g_all[i], ",")) 
        gs = string.(split(gs_all[i], ","))

        println(line)
        println( "BOOTSTRAPPING: ", q_all[i] ) 
        println(line)
        
        cpt, tt, f, adjM, adjMB, probout = bne("BN.data", "BN.header", algo="sm",  scoring_method="BIC", plot="net", DAG=false,  f=f, fs=fs, g=g, gs=gs, relrisk=true, rr_bootstrap=rr_bootstrap, bootout=bootout)    
        
    end

    println("Bootstrap runs finished.")
    println("Results written to file: $bootout")

end
