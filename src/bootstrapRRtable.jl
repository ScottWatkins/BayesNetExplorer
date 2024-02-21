"""

    bootstrapRRtable(RRtable, method="network")

Bootstrap the absolute and relative risk ratio estimates from a network or simple conditional query. Input is the RRtable from RRcalculator.

For network bootstrapping, the BN.data and BN.header files must be available in the directory. For cpq bootstrapping, the original data file must be provided. The output network and simple query tables can be merged and exact test p-values corrected for multiple comparisons.

Options
-------
        method          bootstrap method                      "network|cpq"
        data            original data, all variables           dataframe
                        (required for cpq method)
        rr_bootstrap    number of bootstraps to perform        100
        outfile         basename for the outfile              "bootstrap"
        confmeth        fit CI95 to t-dist or empirical       "empirical|t-dist"
        algo            structure learning method             "sm"
        scoring_method  network scoring method                "BIC"
        plot            plot method                           "net|mb"
        DAG             plot DAG                               true
        data            input data (required for cpq)          ""
        mincounts       min target and cond counts for cpq     [0,1]
        digits          decimal places for cpq                 4
        merge           merge network & cpq outfiles           false
        fdr             set fdr for BH p-value correction      0.05
                        (applied in merge only)

Examples
-------
1. Network bootstrap

  dfnb = bootstrapRRtable(dfr, method="network", rr_bootstrap=500, confmeth="empirical", plot="net",  outfile="mynetboots")

2. Conditional probability bootstrap

  dfpb = bootstrapRRtable(dfr, method="cpq", rr_bootstrap=1000, confmeth="empirical", plot="net",  outfile="mycpqboots", data="myrawdata.csv")

3. Combine results from 1 and 2 correcting for multiple tests

  bootstrapRRtable(dfr, merge=true) 

"""
function bootstrapRRtable(RRtable::Union{String,DataFrame}; method="network", algo="sm", scoring_method="BIC", plot::String="net", DAG::Bool=true, rr_bootstrap::Int=100, outfile::String="bootstrap", merge::Bool=false, confmeth::String="empirical", data::Union{DataFrame,String}="", mincounts::Array{Int,1}=[0,1], digits::Int64=4, fdr=0.05)
    
    if merge == true
        if isfile("$outfile.network.csv") && isfile("$outfile.cpq.csv")
            df1 = CSV.read("$outfile.network.csv", DataFrame);
            df2 = CSV.read("$outfile.cpq.csv", DataFrame);
            dfm = innerjoin(df1, df2, on=[:Query, :Target, Symbol(names(df1)[3]), Symbol(names(df1)[4]), :AbsRiskRatio, :RelRiskRatio, :NonPropProb, :Count, :CondCount, :CondVariables, :CondStateNames, :CondStates,])

            println("Added BH multiple-testing p-value correction with fdr of $fdr...")

            bhc_bin_pval = benjhoc(dfm.Binomial_pval, verbose=false, fdr=fdr)
            bhc_ARR_f_pval = benjhoc(dfm.ARR_Fisher_2t_pval, verbose=false, fdr=fdr)
            bhc_ARR_right_pval = benjhoc(dfm.ARR_Fisher_1t_right, verbose=false, fdr=fdr)
            bhc_RRR_f_pval = benjhoc(dfm.RRR_Fisher_2t_pval, verbose=false, fdr=fdr)
            bhc_RRR_right_pval = benjhoc(dfm.RRR_Fisher_1t_right, verbose=false, fdr=fdr)
            
            insertcols!(dfm, :BH_Binomial_pval => bhc_bin_pval)
            insertcols!(dfm, :BH_ARR_Fisher_2t_pval => bhc_ARR_f_pval)
            insertcols!(dfm, :BH_ARR_Fisher_1t_right => bhc_ARR_right_pval)
            insertcols!(dfm, :BH_RRR_Fisher_2t_pval => bhc_RRR_f_pval)
            insertcols!(dfm, :BH_RRR_Fisher_1t_right => bhc_RRR_right_pval)
            
            CSV.write("$outfile.merged.csv", dfm)
            println("Merged file written to $outfile.merged.csv")

            return dfm

        else
            error("\n\nCould not find the reqired files $outfile.network.csv\nand $outfile.cpq.csv.\n\nPlease run the network and cpq bootstrap methods individually\nto create the two files outfile.network.csv and outifle.cpq.csv,\nthen use the merge=true argument to combine the two files.\n\n")
        end    

    end

    if method == "network"

        if isfile("BN.data") && isfile("BN.header")

            netnodes = readline("BN.header")
            println("Found BN.data and BN.header...\nNetworks for bootstraping will contain: ", netnodes)
            outfile = basename(outfile) * ".network.csv"
            isfile(outfile) ? rm(outfile) : nothing

            OUT = open(outfile, "w")
            println(OUT, "Target\tConditionals\tRR_Conditionals\tP(T)\tP(T|C)\tNet_ARR\tNet_ARRest\tNet_ARRlower\tNet_ARRupper\tNet_RRR\tNet_RRRest\tNet_RRRlower\tNet_RRRupper")
            close(OUT)

        else
            error("BN.data and BN.header files needed. Please run format_file() to recreate these files.")
        end
        
    elseif method == "cpq"
        
        outfile = basename(outfile) * ".cpq.csv"
        isfile(outfile) ? rm(outfile) : nothing
       
        if typeof(data) == DataFrame
            data = data
        elseif isfile(data) && length(data) > 1 
            data = CSV.read(data, DataFrame, types=String);
            println("Found disk file $data...\n")
        else
            error("\n\nOriginal data file or dataframe needed for cpq bootstrapping.\n\n")
        end
    end
            
    line="-"^70
    
    if typeof(RRtable) == DataFrame
        dfr = RRtable
    elseif isfile(RRtable)
        dfr = CSV.read(RRtable, DataFrame, types=String);
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
    qlist = []

    for i in eachindex(f_all)
        
        f = string(f_all[i])
        fs = string(fs)
        g = string.(split(g_all[i], ",")) 
        gs = string.(split(gs_all[i], ","))

        gstr = [ g[i] * "=" * gs[i] for i in eachindex(g) ] 
        gstr = join(gstr, ",")

        query = "P($f=$fs|" * gstr * ")"
        push!(qlist, query)
        
        println(line)
        println( "Bootstrapping: ", q_all[i] ) 
        println(line)
        
        
        if method == "network"
            
            cpt, tt, f, adjM, adjMB, probout = bne(query; algo=algo, scoring_method=scoring_method, plot=plot, DAG=false, relrisk=true, rr_bootstrap=rr_bootstrap, bootout=outfile, confmeth=confmeth)    

        elseif method == "cpq"

            cpq(data, query, bootstraps=rr_bootstrap, confmeth=confmeth, mincounts=mincounts, outfile=outfile, binomial=true, fisher=true, digits=digits)
            
        else
            error("\n\nPlease select a method for bootstrapping (network|cpq)\n\n")
        end
        
    end   #end bootstraps

    println("Bootstrap runs finished.")

    if method == "network"
        
        println("Combining the RRtable and the bootstrap estimates...")
        dfb = CSV.read(outfile, DataFrame, delim="\t");
        dfm = hcat(dfr, dfb[:, [7,8,9,11,12,13]])
        insertcols!(dfm, 1, :Query => qlist)
        CSV.write(outfile, dfm)
        println("Final results written to $outfile.")

    elseif method == "cpq"

        println("Combining the RRtable and the bootstrap estimates...")
        run(`bash -c "cat header.txt $outfile > j.tmp; rm header.txt; mv j.tmp $outfile"`)
        dfc = CSV.read(outfile, DataFrame);
        dfm = hcat(dfr, dfc[:, 4:end])
        insertcols!(dfm, 1, :Query => qlist)
        CSV.write(outfile, dfm)
        println("Final results written to $outfile.")
    end
        
    return dfm
    
end
