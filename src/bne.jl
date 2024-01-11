"""
    dfp, vf, tf, adjM, adjMb, p = bne("query"; options...);

Analyze and explore data using a Bayesian network.

First, create input file with format_file(). Input data can be binary, discrete, or discretized for continuous data [1,2, ..., N].

Return values are: 1) a dataframe of query probabilities, 2) a dataframe of all variable frequencies, 3) target frequencies, 4) adjacency matrix, 5) adjacency matrix for the target Markov blanket, and 6) a probabilities array.

    Network learning methods:

      algo="sm"
        sm   (Silander-Myllymaki exact search)
        hc   (Hill Climbing method)
        mmhc (Max-Min Hill-Climbing heuristic, default)
        sem  (Structural Expectation-Maximization, EM based)

    Scoring methods:

      scoring_method="BIC"
        BDeu (Bayesian-Dirichlet equivalent uniform, default)
        AIC  (Akaike Information Criterion)
        BIC  (Bayesian Information Criterion)
    
    Probability options:

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

    Misc options:

        minfreq             min variable state frequency cutoff [0.0]
        verbose             verbose output                      [false]
        confmeth            CI95 method              [t-dist|empirical]
        rrrdenom            specify target states for the       []
                            relative risk ratio denominator.
                            Default: all query states negated.
        type                propagated probability query type   [conditional]
                            (also joint, marginal; see below)
Examples:

    Step 1. Use format_file() to create the input files.

    Using 2 to 15 non-collinear conditionally dependent variables with minimum frequencies of 3-5 percent (recommended).

    df, ids, traits = format_file("data.csv")

    Step 2. Run the BayesNetExplorer

    1. A simple query using an exact net to estimate the probability of hypertention (HTN) given two evidence variables and states:

    bne( "P(HTN=Yes|Sex=Male,BMI=Low)" )

    2. Query a target with two specified evidence variables using a hill-climb algorithm with AIC scoring. Return probabilies of the target for all combinations of the conditional variable and their states.

    dfp, vf, tf, adjM, adjMb, p = bne("P(A=a|B=b,C=c)", algo="hc", scoring_method="AIC", iterate="nodes" )

    3. Query a target for all variables (nodes) in all states in the network.

     dfp, vf, tf, adjM, adjMb, p = bne("P(HTN=Yes)", algo="hc", scoring_method="BDeu", iterate="all" )

For networks larger than 15 nodes, please use non-exact algorithms (e.g., hc) and set nolimit=true. Network creation time increases rapidly.

Important: Multi-target queries are possible, but input order is imporant! P(A=a,B=b|C=c) is different than P(B=b,A=a|C=c). The conditional propogated query evaluates each variable state left to right. That is, P(B=b,A=a|C=c) will condition on C=c. Then, A=a is evaluted. The probability of B=b in those sample with A=a is returned. Most users will actually want the joint probability query, selected by the kwarg type="joint", which returns the joint probability of A=a and B=b, P(A=a,B=b), conditioned on C.

In all cases, input should be informed by the DAG and the direction of impact of the test variables on the target. Be sure to look at the DAG and the Markov blanket of the targets.
"""
function bne(query::String, data::String="BN.data", header::String="BN.header"; algo="sm", scoring_method="BIC",  bootstrap=0, impute::Bool=false, iterate::String="", minfreq::Float64=0.00, verbose::Bool=false, relrisk::Bool=false, rr_bootstrap::Int=100, bootstrap_method::String="resample", bootout="",  DAG::Bool=false, plot::String="net", boot_plot::Bool=false, nolimit::Bool=false, confmeth::String="t-dist", plotmethod::Symbol=:stress, rrrdenom::Array=[], type="conditional")

    #query parser now replaces f,fs,g,gs manual input
    atargs, aconds, f, fs, g, gs = queryparser(query)

    dfcleaned=DataFrame()

    if isfile("BN.data") && isfile("BN.header")
        println("Found input files BN.data, BN.header...")
    else
        error("\n\nCan't find input files!\nPlease use format_file to format your data!\n\n")
    end
    
    if length(iterate) > 0
        
        if !occursin(r"(^all$)|(^nodes$)", iterate)
            error("\n\nIterate options are \"nodes\" or \"all\". Input was $iterate\n\n")
        end

        if length(rrrdenom) > 0 || typeof(f) == Array
            error("\n\nIteration mode: Only single feature allowed; rrrdenom not yet implemented for iteration\n\n")
        end
        
        if relrisk == true
            error("\n\nIterate mode is not compatible with relrisk.\n\n")
        end
    end

    if occursin(r"\b(sm|hc|mmhc|sem)\b", algo )
        println("Building Bayes net using the $algo algorithm and the $scoring_method scoring method ...")
    else
        error("\n\nERROR: $algo is not an implemented algorithm.\n\n")
    end

    if confmeth == "t-dist" || confmeth == "empirical"
        println("Using $confmeth distribution for CI95...")
    else
        error("\n\nInput confmeth=\"$confmeth\" not recognized.\nOnly t-dist or empirical methods implemented.\n\n")
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

    print("Checking input features...")
    if typeof(f) == Vector{String}
        if length(f) != length(fs)
            error("Each target must have a state, saw $f and $fs")
        end
        for i in f
            if in(i, headdat) == false
                error("\n\nTarget $i not found.\n\n")
            end
        end
    elseif typeof(f) == String
        if in(f, headdat) == false
            error("\n\nTarget $f not found.\n\n")
        end
    end

    println("OK")
    println("Target feature: $f")
    print("Getting variable frequencies...\n")

    # Full tables only; input whole df should be made with format_file.jl
    dfcleaned, df_freq = clean_dfvars_by_frequency("$data", minfreq; nolabels=true, delim=" ", header=headdat)
    global dfcleaned
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
                println("Conditional variable: ", g[i], " state: ", gs[i], " ... OK")
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

    boots = rcopy(R"b <- has.boots(dataset)")
    bootcount = rcopy(Int.(R"bc <- num.boots(dataset)"))
   
    println("Creating the network...")

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
    
    smooth = 1/size(dfcleaned, 1)         #set smooth to reflect dataset size

    @suppress R"df = read.csv(file = $data, sep = ' ')"      # get data table
    @suppress R"colnames(df) <- vars"                        # label data table
    @suppress R"ecpt <- extractCPT(df, dadjm, smooth=$smooth)"  # gRain extract list of CPTs and smooth the NaNs  
    @suppress R"pt <- compileCPT(ecpt)"                      # compile with gRain
    @suppress R"pnet <- grain(pt, smooth=$smooth)"           # make net compiled, smooth
    @suppress R"pnet <- propagate(pnet)"                     # make net propagated
                                                             # pnet goes into ConProb function

    println("Network compilation is...", rcopy(R"isCompiled(pnet)"))
    println("Network propagation is...", rcopy(R"isPropagated(pnet)"))

    CPT = rcopy(R"cpt <- cpts(net)");     #bnstruct conditional prob tables

    # Rplot of DAG
    # plotout=rcopy(R"plotobj <- plot(net, method='qgraph', label.scale.equal=T, node.width = 1.6, plot.wpdag=F)")

    dagcolnames = rcopy(R"colnames(df) <- vars")
    adjM = Int.(adjM)

    if typeof(f) == String     # handle multi-target input
        ff = [f]
    else
        ff = f
    end
    
    mbM, mbnames = get_markov_blanket(adjM, dagcolnames, ff)
    
    if plot == "mb" && sum(mbM) < 1        
        error("\n\nA Markov blanket for $f was not found current network variables.\nPlease set plot=\"net\" or try using other variables.\n\n")
    elseif plot == "net" && sum(adjM) == 0
        error("\n\nNo network connections could be inferred for the input variables.\nPlease try other variables.\n\n")
    end
    
    if plot == "mb"
        mb_plt = plot_network(mbM, headerfile="BN.header", fnode=ff, gnodes=g, DAG=DAG, method=plotmethod )
        println("Displaying Markov blanket plot of $ff ...")
        println("Nodes for minimal Markov blanket...")

        for i in eachindex(ff)
            mbl = join(mbnames[i], ",")
            mbn = ff[i]
            printstyled("$mbn: $mbl\n", color=:cyan)
        end

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
    #println("List of most probable state for each node...")
    #for (k,v) in nv
    #    print(k, ":", v, ";  ")
    #end
    
    df_map = CSV.read("recoded.map", DataFrame, delim=",")
    sort!(df_map, [:feature, :numstate])
    df_target = DataFrame()

    if typeof(f) == Vector{String}
        for i in f
            df_ti =  df_freq[df_freq.feature .== i, :]
            df_target = vcat(df_target, df_ti)    
        end
    else 
        df_target = df_freq[df_freq.feature .== f, :] 
    end

    dftf = df_target.percent .* 0.01 #to frequency

    if length(dftf) == 0
        error("Requested target not found. Check input (f=\"target\")!")
    end

    df_j = leftjoin(df_freq, df_map, on=[:feature, :numstate])
    sort!(df_j, [1,2])
    
    printstyled("\n\nList of variable names, state and frequencies...", color=:green)
    println(df_j, "\n")
    
    proball = Dict()
    mc = 1
    
    println(line)
    println("Calculating probabilities...")

    if iterate == "all"

        if length(g) > 0 || length(gs) > 0
            printstyled("INFO: Iterating all nodes in the network.\nIgnoring user-specified condititional states.\n", color=:yellow)
            println(line)
        end
        
        if length(vars) > 12
            error("Final variable count is >12.\nConsider only the key variables to iterate over\nto complete in reasonable time!")
        end
        
        if length(f) == 0
            error("You must set the target feature for iteration (f=\"target\").")
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

            probout, jpout = ConProb(; f=f, fs=fs, g=g, gs=gs, vars=vars, verbose=verbose, rr_bootstrap=0, type=type, query=query)
            proball[mc] = jpout
            mc += 1
            probout2 = probout

        end

    elseif iterate == "nodes"

        if length(g) < 2
            error("\n\nPlease specify at least 2 conditional variables to test for iteration.\n\n")
        end

        println(line)
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
            
            probout, jpout = ConProb(; f=f, fs=fs, g=g, gs=gs, vars=vars, verbose=verbose, rr_bootstrap = 0, type=type, query=query)
            proball[mc] = jpout
            mc += 1
            probout2 = probout 
        end

    else


        println("Performing probability query...\nListing target $type probability estimates...")

        if length(rrrdenom) > 0
            println("Using user-specified relative risk ratio denominator states: $rrrdenom")
            gso = rrrdenom
        else
            gso = replace(gs, "1" => "2", "2" => "1")
        end

        fsi = parse.(Int64, fs)
        gsi = parse.(Int64, gso)

        function mtarg_freq(f, fs)  #Get multi-targ freq

            qin = ""
            for i in eachindex(f)
                q = "dfcleaned." * string(f[i]) * " .== " * string(fs[i]) * ","
                qin = qin * q
            end
            qin = qin[1:end-1]
            qin = "dfcleaned[ .&( " * qin * "), :]"
            qin = Meta.parse(qin)
            return qin
        end
        
        if typeof(f) == Vector{String}
            qin = mtarg_freq(f, fs)
            dfqin = eval(qin)
            cqin = size(dfqin,1)
            mtfreq = cqin/ size(dfcleaned,1)
            tpf = mtfreq    # obs. baseline freq multi-target
            println("Observed $cqin samples with ", join(atargs, ", "), "\nBaseline frequency: ",  round(tpf, digits=6), "\n", line)
        else
            tpf = dftf[fsi] # obs. baseline freq single target
            println("Observed baseline frequency of ", join(atargs, ", "), ": ",  round(tpf, digits=6), "\n", line)
        end

        probout, jpout = ConProb(; f=f, fs=fs, g=g, gs=gs, vars=vars, verbose=true, rr_bootstrap=0, type=type, query=query)
        
        probout_o, jpout_o = ConProb(; f=f, fs=fs, g=g, gs=gso, vars=vars, verbose=false, rr_bootstrap=0, type=type, query=query);        

        if type == "conditional" || type == "joint"
            rrr = probout[CartesianIndex(Tuple(fsi))]
        
            rro = probout_o[CartesianIndex(Tuple(fsi))]     # rr opp prob, idx still fsi!

            RR = round(rrr / rro, digits=6)                 # This is the init RR est

            absrisk = round(rrr/tpf, digits=6)

        end  
            
        
        if relrisk == true

            type == "marginal" ? error("\n\nBootstrapped relative risk not implemented for $type proabilities.\n\n") : ""

            println(line)
            printstyled("Calculating relative and absolute risks with CI95 ($rr_bootstrap bootstraps) ...\n", color=:green)
            
            CI95, PRdist, RRest = bne_bootstrap(data, header; impute=false, algo=algo, scoring_method=scoring_method, f=f, fs=fs, g=g, gs=gs, bootstrap=0, iterate="", minfreq=0.00, verbose=false, rr_bootstrap=rr_bootstrap, bootstrap_method=bootstrap_method, boot_plot=boot_plot, confmeth=confmeth, rrrdenom=rrrdenom, type=type, query=query)

            ARdist  = sort(PRdist ./ tpf)

            if (confmeth == "t-dist")

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

            if typeof(f) == Vector{String}
                fcond = join([f[i] * "=" * string(fs[i]) for i in eachindex(f)], "," )
                oconds = join([g[i] * "=" * string(gso[i]) for i in eachindex(g)], ", ")
            else
                fcond = f * "=" * fs
                oconds = join(g .* "=" .* gso, ", ")
            end
            
            println("\n", line, "\n")
            rrr = round(rrr, digits=6)
            rro = round(rro, digits=6)

            println("Absolute Risk Estimates:\n")
            println("    P($fcond|$pconds) = $rrr")
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
            
            if (absrisk < ARlower) || (absrisk > ARupper) || (RR < CI95[1]) || (RR > CI95[2])
                printstyled("\nWARNING:\nOne or more network risk estimates are outside the risk distribution.\n", color=:yellow)
                printstyled("This can be due to unstable network branches, small sample size, etc.\nPlease condsider node selection, sample size, and the probability query.\n", color=:yellow)
            end
            
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
        probout2 = probout
        probout = join([f, probout, g, gs], "|")
        proball[mc] = probout
        
    end

    if type == "marginal"
        return "NA", "NA", "NA", "NA", "NA", "NA"
    end
    

    if typeof(f) == Vector{String}  #if complex target don't write file for more analysis

        fsn = Vector(parse.(Int64,fs))
        pF = round(probout2[CartesianIndex(Tuple(fsn))], digits=4)
        if type == "joint"  # clarify conditional multi-targets results  
            type = type * " target"
        end
        printstyled("Network $type probability estimate of\n$query = $pF\n", color=:green)
        println(line)
        
        
        h = [:Targets, :ProbabilityMatrix, :ConditionalVariable, :ConditionalStates]
        dfp = DataFrame([name => [] for name in h])
        mt =  string.(split(replace(probout, "\"" => ""), "|"))
        push!(dfp,mt)
        insertcols!(dfp, 2, :States => string(join(fs, ",")))
        insertcols!(dfp, 3, query => pF)
        
        return dfp, df_j, dftf, adjM, mbM, probout #cpdf, df of used variables, target state freq
        
    else

        fsn = parse(Int64,fs)
        pF = round(probout2[CartesianIndex(Tuple(fsn))], digits=4)
        printstyled("Network $type estimate of $query: ", pF, "\n", color=:teal)
        println(line)
        
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

        return dfp, df_j, dftf, adjM, mbM, probout #cpdf, df of used variables, target state freq

    end
end
