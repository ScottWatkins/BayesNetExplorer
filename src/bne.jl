"""
    cpt, tt, tf, adjM, adjMb, probout = bne("BN.data", "BN.header"; algo="sm", scoring_method="bic", f="", fs=0,  g=[], gs=[])

    Construct a Bayesian network. Create the input BN.data and BN.header files with format_file(). Data can be binary, discrete, or discretized for continuous data [1,2, ..., N]. Sample names are always in column1, and data should be type consistent, that is, all strings or all bool.

Return values include  1) dataframe of conditional probabilities, 2) a dataframe of all variable frequencies, 3) target variable frequencies 4) adjacency matrix, 5) adjacency matrix for Markov blanket, and 6) a probabilities array. Additional options can be specified using keywords.

    Network learning methods:

        sm   (Silander-Myllymaki exact search)
        hc   (Hill Climbing method)
        mmhc (Max-Min Hill-Climbing heuristic, default)
        sem  (Structural Expectation-Maximization, EM based)

    Scoring methods:

        BDeu (Bayesian-Dirichlet equivalent uniform, default)
        AIC  (Akaike Information Criterion)
        BIC  (Bayesian Information Criterion)
    
    Probability options with input example:

        f                Feature/Target variable (required)     "F1"
        g                Given/Evidence variables               ["F2", "F3"]
        gs               Given/Evidence  states                 ["2",   "1"]
                         (Must be numeric string)
        iterate          Iterate states                         [nodes|all]
        relrisk          calc. rel risk ratio                   [false]

                     
    Bootstrapping options

        rr_bootstrap        number bootstrap interations        [100]
        bootstrap_method    resample or delete half             [resample]
        bootout             write and append results to file    ""
        boot_plot           show plot for each bootstrap        false

    Plot options

        DAG                 plot is a directed acyclic graph    [false]
        plot                network or markov blanket "mb"      ["net"]
        plotmethod          network plot style                  [:stress]

    Other options:

        minfreq             min variable state frequency cutoff [0.0]
        verbose             verbose output                      [false]
        confmeth            CI95 calculation method             [t-dist|normal]
        rrrdenom            specify target states for the       []
                            relative risk ratio denominator.
                            Default: all query states negated.
    Examples:

    First, use format_file() to create the input files for the network. For exact networks, use 5 to 15 non-collinear conditionally dependent variables with minimum frequencies of five percent (recommended). For larger networks, use non-exact algorithms and set nolimit=true.

    1. A simple query for a feature and two given evidence variables and states:
    bne("BN.data", "BN.header", algo="sm", scoring_method="BIC", f="A", g=["B", "C" ], gs=["2","1"] )

    2. Query feature for 3 specified evidence variables. Examine all states for all evidence nodes.
    cpt, tf, f, adjM, adjMb, probout = bne("BN.data", "BN.header", algo="hc", scoring_method="AIC", f="A", g=["B", "C", "D" ], iterate="nodes" )

    3. Query a feature for variables (nodes) in all states in the network.
    cpt, tf, f, adjM, adjMb, probout = bne("BN.data", "BN.header", algo="hc", scoring_method="BDeu", f="A", iterate="all" )


"""
function bne(data, header; algo="sm", scoring_method="BIC", f::String="", fs=0, g::Array=[], gs::Array=[], bootstrap=0, impute::Bool=false, iterate::String="", minfreq::Float64=0.00, verbose::Bool=false, relrisk::Bool=false, rr_bootstrap::Int=100, bootstrap_method::String="resample", bootout="",  DAG::Bool=false, plot::String="net", boot_plot::Bool=false, nolimit::Bool=false, confmeth::String="normal", plotmethod::Symbol=:stress, rrrdenom::Array=[])
    
    if length(iterate) > 0
        
        if !occursin(r"(^all$)|(^nodes$)", iterate)
            error("\n\nIterate options are \"nodes\" or \"all\". Input was $iterate\n\n")
        end
        
        if relrisk == true
        error("\n\nIterate option is not compatible with relrisk.\n\n")
        end

    end

    if occursin(r"\b(sm|hc|mmhc|sem)\b", algo )
        println("Building Bayes net using the $algo algorithm and the $scoring_method scoring method ...")
    else
        error("\n\nERROR: $algo is not an implemented algorithm.\n\n")
    end

    line = "-"^75   
    
    @suppress R"library(bnstruct)"
    @suppress R"library(qgraph)"
    @suppress R"library(gRain)"
    
    println("Input data file: $data")
    println("Input header file: $header")
    vardat = readlines("$header")
    headdat = replace(readline("$header"), "\"" => "")
    headdat = String.(split(headdat, r"\s+"))

    if length(bootout) > 0 && !isfile("bootheader.txt")   #header for bootstrap file
        OUT = open("bootheader.txt", "w")
        println(OUT, "## Network nodes: $headdat ")
        println(OUT, "Target\tConditionals\tCondDenom\tP(T)\tP(T|C)\tARR\tARRest\tARRlower\tARRupper\tRRR\tRRRest\tRRRlower\tRRRupper")
        close(OUT)
    end
    
    limit = 16
    if nolimit == true
        limit = 10000
    else
        if length(headdat) > limit
            error("\n\nBNE: Input with >16 variables will have very long compute times.\nSelect up to 16 variables for exact network analysis.\nFor example, format_file(\"filename.csv\", datacols=[1,2,5,6,7]).\n\nLarge approximate networks can be made with using\ngreedy algorithms (e.g., algo=\"hc\") and setting nolimit=true.\nGreedy networks may differ from the more accurate exact networks.\n\n")
        end
    end

    numstates = parse.(Int64, split(vardat[2]))
    types = split(vardat[3])

    if in(f, headdat) == false
        error("\n\nTarget $f not found.\n\n")
    end

    println("Target feature: $f")
    print("Getting variable frequencies...\n")

    # Full tables only; cleaned df should be made with format_file.jl
    dfcleaned, df_freq = clean_dfvars_by_frequency("$data", minfreq; nolabels=true, delim=" ", header=headdat)

    #check input and features
    if length(iterate) > 0
        
    else
        println("Checking conditional variables and states...")

        if sum(occursin.(r" ", gs)) > 0
            error("\nSpaces are not allowed in the gs states: gs = $gs\n\n")
        end
        
        if length(g) != length(gs)
            error("\nPlease provide one state for each conditional variable.\ng:  $g\ngs: $gs\n")
        end

        for i in eachindex(g)
            if in(parse(Int64, gs[i]), dfcleaned[:, Symbol(g[i])])
                println("Conditional variable: ", g[i], " state: ", gs[i], " ... OK.")
            else
                error("\n\nConditional variable: ", g[i], " state: ", gs[i], " ... not found. Check input.\n\n")
            end
        end
    end
        
    print("Processing network data...\n")
    
    R"dataset <- BNDataset( $data, $header, starts.from = 1)"  #must use 1-based variables

    vars = rcopy(R"vars <- variables(dataset)")
    exv = length(vars) - size(dfcleaned, 2) 

    if exv > 0
        error("Found $exv variables that do not meet a minfreq of $minfreq.\nPlease rerun format_files() to filter data.")
    end
    
    if impute == true
        error("\n\nImputation is now implemented in the impute_dataframe() function\n")
        #R"dataset <- impute(dataset)"
    end

    # use @suppress to suppress screen output by R 

    if impute == true && bootstrap > 0
        println("Using imputed data and network bootstrapping $bootstrap times...")
        @suppress R"dataset <- bootstrap(dataset, num.boots = $bootstrap, imputation = $impute)"
        @suppress R"net <- learn.network(dataset, algo = $algo, scoring.func = $scoring_method, use.imputed.data = $impute)"
    elseif impute == false && bootstrap > 0
        println("Using fully known data and network bootstrapping $bootstrap times...")
        @suppress R"dataset <- bootstrap(dataset, num.boots = $bootstrap)"
        @suppress R"net <- learn.network(dataset, algo = $algo, scoring.func = $scoring_method)"
    elseif impute == true && bootstrap == 0
        @suppress println("Using imputed data and without network bootstrapping...")
        @suppress R"net <- learn.network(dataset, algo = $algo, scoring.func = $scoring_method, use.imputed.data = $impute)"
    elseif impute == false && bootstrap == 0
        @suppress println("Using fully known data without network bootstrapping...")
        @suppress R"net <- learn.network(dataset, algo = $algo, scoring.func = $scoring_method)"
    end

    println("Network information:")
    
    boots = rcopy(R"b <- has.boots(dataset)")
    bootcount = rcopy(Int.(R"bc <- num.boots(dataset)"))
   
    println("Creating the network ...")

    if boots == true
        @suppress R"wpnet  <- wpdag.from.dag(net)"  #wp net object
        @suppress R"wpadjm <- wpnet@wpdag"          #bootstrapped adj. matrix
        adjM = rcopy(R"dadjm  <- wpnet@dag")            #non-bootstraped adj. matrix        
        @suppress R"rownames(wpadjm) <- vars; colnames(wpadjm) <- vars" #label boot matrix
        @suppress R"rownames(dadjm) <- vars; colnames(dadjm) <- vars"   #label non-boot matrix
    elseif boots == false
        adjM = rcopy(R"dadjm <- net@dag")                #non-bootstrapped DAG adj matrix
        @suppress R"rownames(dadjm) <- vars; colnames(dadjm) <- vars" #label DAG matrix
    end

    @suppress R"df = read.csv(file = $data, sep = ' ')"      # get data table
    @suppress R"colnames(df) <- vars"                        # label data table
    @suppress R"gn <- as(dadjm, 'graphNEL')"                 # DAG matrix
    @suppress R"ecpt <- extractCPT(df, gn, smooth = 0.01)"   # extract CPT from data, smooth NaN  
    @suppress R"pt <- compileCPT(ecpt)"                      # 
    @suppress R"pnet <- grain(pt, propagate=T)"              # make cpt net, propagated, copy
                                                             # pnet goes into CondProb function
    CPT = rcopy(R"cpt <- cpts(net)");     #bnstruct conditional prob tables

#   plotout = rcopy(R"plotobj <- plot(net, method='qgraph', label.scale.equal=T, node.width = 1.6, plot.wpdag=F)")
    dagnames = rcopy(R"colnames(df) <- vars")

    adjM = Int.(adjM)

    mbM, mbnames = get_markov_blanket(adjM, dagnames, f)

    if plot == "mb" && sum(mbM) < 2
        error("\nMarkov blanket for $f not found.\nTry using more variables or setting plot=net.\n")
    end
    
    if plot == "mb"
        mb_plt = plot_network(mbM, headerfile="BN.header", fnode=f, gnodes=g, DAG=DAG, method=plotmethod )
        println("Displaying Markov blanket plot of $f ...")
        println("Nodes for minimal Markov blanket for $f...")
        printstyled("$mbnames\n", color=:cyan)
        display(mb_plt)
        #savefig(mb_plt, "$f.markov.svg")
    elseif plot == "net"
        plt = plot_network(adjM, headerfile="BN.header", fnode=f, gnodes=g, DAG=DAG, method=plotmethod ) #DAG in, converted if UG
        display(plt)
    else
        println("Plot must be either net (entire net) or mb (Markov blanket of target) or skip.")
    end
    
    nodevals = rcopy(Int.(R"get.most.probable.values(net, prev.values = NULL)")) 

    nv = Dict(zip(vars, (nodevals)))

    #println("Most common state for each node in the DAG:")
    #for (k,v) in nv
    #    println(k, "\t", v)
    #end

    df_map = CSV.read("recoded.map", DataFrame, delim=",")
    sort!(df_map, [:feature, :numstate])
    
    df_target = df_freq[df_freq.feature .== f, :] 
    dftf = df_target.percent .* 0.01 #to frequency

    if length(dftf) == 0
        error("Requested target not found. Check input (f=\"target\")!")
    end

    df_j = leftjoin(df_freq, df_map, on=[:feature, :numstate])
    sort!(df_j, [1,2])
    
    printstyled("\n\nList of states, names, and frequencies for all variables...", color=:green)
    println(df_j, "\n")
    
    proball = Dict()
    mc = 1

    println("Calculating probabilities...")
    
    if iterate == "all"

        if length(rrrdenom) > 0
            error("\n\nrrrdenom not yet implemented for iteration\n\n")
        end

        if length(g) > 0 || length(gs) > 0
            printstyled("INFO: Iterating all nodes in the network.\nIgnoring user-specified condititional states.\n", color=:yellow)
        end
        
        if length(vars) > 12
            error("Final variable count is >12.\nConsider only the key variables to iterate over\nto complete in reasonable time!")
        end
        
        if length(f) == 0
            error("You must set the target feature for interation (f=\"target\").")
        end
        
        ivarsall_set = create_powerset(df_freq, f)
        popfirst!(ivarsall_set) #remove empty set
        
        for i in eachindex(ivarsall_set)
            
            #Recompose variables and states
            g, gs = recompose_powerset(ivarsall_set[i]);
            
            #Remove powerset combinations with duplicates
            if length(unique(g)) < length(g)  
                continue
            end

            probout, jpout = ConProb(; f=f, fs=fs, g=g, gs=gs, vars=vars, verbose=verbose, rr_bootstrap=0)
            proball[mc] = jpout
            mc += 1
        end

    elseif iterate == "nodes"

        if length(rrrdenom) > 0
            error("\n\nrrrdenom not yet implemented for iteration\n\n")
        end
        
        println("Calculating all probabilities selected nodes: $g")
        
        df_j = df_j[findall(in(g), df_j.feature), :] 
        
        ivarsall_set = create_powerset(df_j, f)
        popfirst!(ivarsall_set)
        
        for i in g
            if !in(i, df_j.feature)
                error("Requested node $i not available")                
            end
        end
        
        for i in eachindex(ivarsall_set)
            
            g, gs = recompose_powerset(ivarsall_set[i])
            
            if length(unique(g)) < length(g)
                continue
            end
            
            probout, jpout = ConProb(; f=f, fs=fs, g=g, gs=gs, vars=vars, verbose=verbose, rr_bootstrap = 0)    
            proball[mc] = jpout
            mc += 1
 
        end

    else

        println(line)
        printstyled("Performing single probability query ...\nListing conditional probabilities for all states of $f ...\n", color=:green)

        if length(rrrdenom) > 0
            println("Using user specified relative risk ratio denominator: $rrrdenom")
            gso = rrrdenom
        else
            gso = replace(gs, "1" => "2", "2" => "1")
        end
        
        probout, jpout = ConProb(; f=f, fs=fs, g=g, gs=gs, vars=vars, verbose=true, rr_bootstrap=0)
        #println("===>", probout, "  ", jpout)

        probout_o, jpout = ConProb(; f=f, fs=fs, g=g, gs=gso, vars=vars, verbose=false, rr_bootstrap=0);        
        
        if relrisk == true

            printstyled("Calculating relative and absolute risks with CI95 ($rr_bootstrap bootstraps) ...\n", color=:green)
            println(line)
            
            if fs == 0 || typeof(fs) == Int
                error("\n\nPlease set the target state to bootstrap in string format (e.g. fs=\"2\").\n")
            elseif in( parse(Int64,fs), dfcleaned[:, Symbol(f)] )
                println("Target: ", f, ", state: ", fs, " ... OK.")
            else
                error("\n\nTarget: ", f, ", state: ", fs, " ... not found. Check input.\n\n")
	    end
            
            fsi = parse(Int,fs)
            rrr = probout[fsi]       #rr target prob
            rro = probout_o[fsi]     #rr opp prob
            
            RR = round(rrr / rro, digits=4)
            
            tpf = dftf[fsi]          #target (baseline) pop freq

            absrisk = round(rrr/tpf, digits=4)

            CI95, PRdist, RRest = bne_bootstrap(data, header; impute=false, algo=algo, scoring_method=scoring_method, f=f, fs=fs, g=g, gs=gs, bootstrap=0, iterate="", minfreq=0.00, verbose=false, rr_bootstrap=rr_bootstrap, bootstrap_method=bootstrap_method, boot_plot=boot_plot, confmeth=confmeth, rrrdenom=rrrdenom)

             
            
            ARdist  = sort(PRdist ./ tpf)

            if (confmeth == "normal")

                ARest = round(median(ARdist), digits=4)
                ARstd = std(ARdist)
                ARlower, ARupper = round.(confint(OneSampleTTest(ARest, ARstd, length(ARdist), 0.05 )), digits=4)

                if ARlower < 0
                    ARlower = 0.0
                end

                elseif (confmeth == "empirical")

                ARest   = round(median(ARdist), digits=4)            
                ARupper = round(ARdist[ Int(floor(length(ARdist) * 0.975)) ], digits=4)
                ARlower = round(ARdist[ Int(ceil(length(ARdist) * 0.025)) ], digits=4)

            end
            
            tpf = round(tpf, digits = 6)
            RRest = round(RRest, digits=4)
            
            pconds = join(g .* "=" .* gs , ", ")
            fcond = f * "=" * fs
            oconds = join(g .* "=" .* gso, ", ")

            println()
            println(line)
            println()

            println("Absolute Risk Estimates:\n")
            println("    P($fcond|$pconds)")
            println("       vs.")
            println("    P($fcond) = $tpf  (baseline)\n")
            printstyled("    Network propagated Absolute Risk Ratio: $absrisk\n", color=:cyan)
            printstyled("    Bootstrap Absolute Risk Distribution:   $ARest  CI95 ($ARlower, $ARupper)\n", color=:cyan)

            println()
            println("Relative Risk Estimates:\n")
            println("    P($fcond|$pconds) = $rrr")
            println("       vs.")
            println("    P($fcond|$oconds) = $rro\n")
            printstyled("    Network propagated Relative Risk Ratio: $RR\n", color=:cyan)
            printstyled("    Bootstrap Relative Risk Distribution:   $RRest  CI95 $CI95\n\n", color=:cyan)
            
#            if (absrisk < ARlower) || (absrisk > ARupper) || (RR < CI95[1]) || (RR > CI95[2])
#                printstyled("\nINFO:\nRisk estimate(s) from input data fall outside the resampled distribution.\n", color=:grey)
#                printstyled("This often means there are few, if any, samples for the requested query\nor that too few bootstraps were performed.\n", color=:yellow)
#            end
            
            println(line)

            if length(bootout) > 0
                tvl = f * "=" * fs
                cvl = ""
                for i in eachindex(g)
                    cvl = cvl * string(g[i]) * "=" * string(gs[i]) * ","
                end

                cvl = cvl[1:end-1]
                gso_p = join([g[i] * "=" * gso[i] for i in eachindex(g)], ",") # printable denom states

                RRlower = CI95[1]
                RRupper = CI95[2]
                OUT = open(bootout, "a")
                println(OUT, "$tvl\t$cvl\t$gso_p\t$tpf\t$rrr\t$absrisk\t$ARest\t$ARlower\t$ARupper\t$RR\t$RRest\t$RRlower\t$RRupper")
                close(OUT)

            end
            
        end

        g = join(g, ",")            #create the proball array if doing a single query
        gs = join(gs, ",")          #for output into the cpt table. Warning g and gs are stringified.
        probout = join([f, probout, g, gs], "|")
        proball[mc] = probout
        
    end
        
    header = ["Target"]     #create header for multiple states of target
    sck = sum(df_map.feature .== f)
    sc = df_map.state[df_map.feature .== f]

    if sck == 0
        error("\n\nFeature $f not found in the current features array.\nIf you changed the BN.header file manually, put the new names in the features array.\n\n")
    end
    
    for i in sc
        ts =  join([ f, "_", i ], "")
        push!(header, ts)

    end

    push!(header,  "ConditionalVariables", "ConditionalStates")

    OUT2 = open("BN.grain.cpts.out", "w") #save bnstruct cpts to disk file
    for i in 1:length(proball)
        v = proball[i]
        v = replace.(v, "[" => "", "]" => "")
        v = replace(v, "," => "|", count = (sck - 1))
        println(OUT2, v)
    end    
    close(OUT2)

    R"rm(list = ls())"    #clear R workspace

    dfp = CSV.read("BN.grain.cpts.out", DataFrame, delim="|", header=header)

    dfp = sort(dfp, [2,3,4], rev=[true,true,true])

    CSV.write("BN_cpts.table", dfp)

    rm("BN.grain.cpts.out")

    return dfp, df_j, dftf, adjM, mbM, probout  #cpdf, df of used variables, target state freq

end
