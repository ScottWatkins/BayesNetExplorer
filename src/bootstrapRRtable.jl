"""

    bootstrapRRtable(RRtable)

Bootstrap the output of the RRcalculator using bne. Previously created BN.data and BN.header files must exist. The variables in the BN disk files define network nodes. 

    Options:
        rr_bootstrap    number of bootstraps to perform        100
        bootout         write results to outfile               "bootstrap.out"
        confmeth        fit CI95 to t-dist or empirical        "t-dist|empirical"
        merge           merge ARR & RRR results to RRtable     true
"""
function bootstrapRRtable(RRtable::Union{String,DataFrame}; rr_bootstrap::Int=100, bootout::String="bootstrap.out", merge::Bool=true, confmeth="t-dist")

    #addprocs(threads)
    if isfile("BN.data") && isfile("BN.header")
        netnodes = readline("BN.header")
        println("Found BN.data and BN.header...\nNetworks for bootstraping will contain: ", netnodes)
    end
    
    if isfile(bootout)
        println("Found existing file called $bootout...")
        println("Removed the existing file...")
        rm(bootout)
    end
    
    OUT = open(bootout, "w")
    println(OUT, "Target\tConditionals\tRR_Conditionals\tP(T)\tP(T|C)\tARR\tARRest\tARRlower\tARRupper\tRRR\tRRRest\tRRRlower\tRRRupper")
    close(OUT)
                
    line="-"^70
    
    if typeof(RRtable) == DataFrame
        dfr = RRtable
    elseif isfile(RRtable)
        dfr = CSV.read(RRtable, DataFrame, types=String)
    end

    q_all = querywriter(dfr)

    f_all = String.(dfr[!,1])

    fs = String(split(names(dfr)[2], r"_|\|" )[2])

    if fs == "Yes"
        println("Setting fs=Yes to fs=2...")
        fs = "2"
    end
    
    g_all = dfr[!,9]
    gs_all = dfr[!,11]

    for i in eachindex(f_all)

        f = string(f_all[i])
        fs = string(fs)
        g = string.(split(g_all[i], ",")) 
        gs = string.(split(gs_all[i], ","))

        println(line)
        println( "BOOTSTRAPPING: ", q_all[i] ) 
        println(line)

        cpt, tt, f, adjM, adjMB, probout = bne("BN.data", "BN.header", algo="sm",  scoring_method="BIC", plot="net", DAG=false,  f=f, fs=fs, g=g, gs=gs, relrisk=true, rr_bootstrap=rr_bootstrap, bootout=bootout, confmeth=confmeth)    

        
    end

    println("Bootstrap runs finished.")
    println("Results written to file: $bootout")

    if merge == true
        println("Merging bootstrap estimates and RRtable...")
        dfb = CSV.read(bootout, DataFrame, delim="\t")
        dfm = hcat(dfr, dfb[:, [7,8,9,11,12,13]])
    end

    return dfm
    
end
