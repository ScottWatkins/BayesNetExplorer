"""
    cpt, tt, tf, adjm = bne("BN.data", "BN.header"; algo="sm", scoring_method="bic", f="", fs=0,  g=[], gs=[], verbose=false, DAG=false)

    Construct a Bayesian network. Inputs are the formatted data and header files (see format()). Data must be binary, discrete or discretized continuous [1,2, ..., N]. Return 1) the conditional probabilities, 2) a dataframe of all input variables, mappings and state frequencies 3) target state frequencies and 4) the graph adjacency matrix.

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

        f            Target feature (required)         "F1"
        g            Given/Evidence features            ["F2", "F3"]
        gs           Given/Evidence  states             ["2",   "1"]
                     (Must be numeric string)
        iterate      Iterate states                     [nodes|all]
        relrisk      calc. rel risk ratio for target    [false]

                     
    Bootstrapping options

        rr_bootstrap        number bootstrap interations        [100]
        bootstrap_method    resample or delete half             [resample]

    Plot options

        DAG    plot is directed acyclic graph    [false]
        plot   network or markov blanket "mb"    ["net"]

    Other options:

        minfreq      min variable state frequency cutoff [0.00]
        verbose      verbose output                      [false]

    Examples:

    1. Simple query for a target and two evidence features and their states:
    bne("BN.data", "BN.header", algo="sm", scoring_method="BIC", f="A", g=["B", "C", ], gs=["2","1"] )

    2. Query feature for 3 specified nodes. All states for all nodes examined.
    cpt, c_df, f, adjM = bne("BN.data", "BN.header", algo="hc", scoring_method="AIC", f="A", g=["B", "C", "D" ], iterate="nodes" )

    3. Query feature for all nodes. All states for all nodes examined.
    cpt, c_df, f, adjM = bne("BN.data", "BN.header", algo="sm", scoring_method="BDeu", f="A", iterate="all" )


"""
function bne(data, header; impute=false, algo="sm", scoring_method="BIC", f::String="", fs=0, g::Array=[], gs::Array=[], bootstrap=0, iterate::String="", minfreq::Float64=0.00, verbose::Bool=false, relrisk::Bool=false, rr_bootstrap::Int=100, bootstrap_method::String="resample", DAG::Bool=false, plot::String="net")

    if sum(occursin.(r" ", gs)) > 0
        error("\nSpaces are not allowed in the gs states: gs = $gs\n\n")
    end

    println("Building bayes net using the $algo algorithm and the $scoring_method scoring method ...")
    
    @suppress R"library(bnstruct)"
    @suppress R"library(qgraph)"
    @suppress R"library(gRain)"
    
    println("Input data file: $data")
    println("Input header file: $header")
    vardat = readlines("$header")
    headdat = replace(readline("$header"), "\"" => "")
    headdat = String.(split(headdat, r"\s+"))

    if length(headdat) > 20
        error("\n\nBNE: input data exceeds the recommended 20 variables. Exiting.\n")
    end
      
    numstates = parse.(Int64, split(vardat[2]))
    types = split(vardat[3])

    if in(f, headdat) == false
        error("\n\nTarget $f not found.\n\n")
    end

    print("Getting variable frequencies...\n")

    # Full tables only; cleaned df should be made with makeBNdataset.jl
    dfcleaned, df_freq = clean_dfvars_by_frequency("$data", minfreq; nolabels=true, delim=" ", header=headdat)

    print("Processing network data...\n")
    
    R"dataset <- BNDataset( $data, $header, starts.from = 1)"  #must use 1-based variables

    vars = rcopy(R"vars <- variables(dataset)")
#    exv = length(vars) - size(dfcleaned, 2) 

#    if exv > 0
#        error("Found $exv variables that do not meet a minfreq of $minfreq.\nPlease rerun format_files() to filter data.")
#    end
    
    if impute == true
        R"dataset <- impute(dataset)"
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
    @suppress R"ecpt <- extractCPT(df, gn, smooth = 0.00001)"  # extract CPT from data, smooth NaN  
    @suppress R"pt <- compileCPT(ecpt)"                      # 
    @suppress R"pnet <- grain(pt, propagate=T)"              # make cpt net, propagated, copy
                                                             # pnet goes into CondProb function
    CPT = rcopy(R"cpt <- cpts(net)");     #bnstruct conditional prob tables

#   plotout = rcopy(R"plotobj <- plot(net, method='qgraph', label.scale.equal=T, node.width = 1.6, plot.wpdag=F)")
    dagnames = rcopy(R"colnames(df) <- vars")
    adjM = Int.(adjM)
    mbM = get_markov_blanket(adjM, dagnames, f)

    if plot == "mb"
        mb_plt = plot_network(mbM, "BN.header", fnode=f, gnodes=g, DAG=DAG )
        println("Displaying markov blanket plot of $f ...")
        display(mb_plt)
        #savefig(mb_plt, "$f.markov.svg")
    elseif plot == "net"
        plt = plot_network(adjM, "BN.header", fnode=f, gnodes=g, DAG=DAG ) #DAG in, converted if UG
        display(plt)
    else
        println("Plot must be either net (entire net) or mb (markov blanket of target) or skip.")
    end
    
    nodevals = rcopy(Int.(R"get.most.probable.values(net, prev.values = NULL)")) 
    nv = Dict(zip(vars, (nodevals)))
    println("Most common state for each node in the DAG:")

    for (k,v) in nv
        println(k, "\t", v)
    end

    df_map = CSV.read("recoded.map", DataFrame, delim=",")
    sort!(df_map, [:feature, :numstate])
    
    df_target = df_freq[df_freq.feature .== f, :] 
    dftf = df_target.percent .* 0.01 #to frequency

    if length(dftf) == 0
        error("Requested target not found. Check input!")
    end

    df_j = leftjoin(df_freq, df_map, on=[:feature, :numstate])
    sort!(df_j, [1,2])
    
    printstyled("\n\nList of states, names, and frequenies for all variables...", color=:green)
    println(df_j, "\n")
    
    proball = Dict()
    mc = 1

    println("Calculating probabilities...")
    
    if iterate == "all"

        if length(g) > 0 || length(gs) > 0
            printstyled("INFO: Iterating all; ignoring user-specified condititional states.\n", color=:yellow)
        end
        
        if length(vars) > 20
            error("Too many variables to iterate all of them!")
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
            
            probout = ConProb(; f=f, g=g, gs=gs, vars=vars, verbose=verbose, rr_bootstrap=0)    
            proball[mc] = probout
            mc += 1
            
        end
        

    elseif iterate == "nodes"

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
            
            probout = ConProb(; f=f, g=g, gs=gs, vars=vars, verbose=verbose, rr_bootstrap = 0)    
            proball[mc] = probout
            mc += 1
            
        end


    else

        println("$('-'^75)")
        printstyled("Performing single probability query ...\nList will include probabilities for all states of $f ...\n", color=:green)
        probout = ConProb(; f=f, g=g, gs=gs, vars=vars, verbose=true, rr_bootstrap=1)        

        gso = replace(gs, "1" => "2", "2" => "1")
        probout_o = ConProb(; f=f, g=g, gs=gso, vars=vars, verbose=false, rr_bootstrap=1);        
        
        if relrisk == true

            if length(unique(gso)) > 2
                error("Relative risk calculations are limited to two-state variable.")
            end
            
            printstyled("Calculating relative and absolute risks with CI95 ($rr_bootstrap bootstraps) ...\n", color=:green)
            println("$('-'^75)")
            
            if fs == 0 || typeof(fs) == Int
                error("\nPlease set indicate the feature-state you wish to bootstrap (\"1\" or \"2\").\n")
            end

            fsi = parse(Int,fs)
            rrr = probout[fsi]       #rr target prob
            rro = probout_o[fsi]     #rr opp prob
            RR = round(rrr / rro, digits=4)
            
            tpf = dftf[fsi]          #target pop freq
            absrisk = round(rrr/tpf, digits=4)
            
            CI95, PRdist = bne_bootstrap(data, header; impute=false, algo=algo, scoring_method=scoring_method, f=f, fs=fs, g=g, gs=gs, bootstrap=0, iterate="", minfreq=0.00, verbose=false, rr_bootstrap=rr_bootstrap, bootstrap_method="resample")
            
            ARdist = sort(PRdist ./ tpf)
            ARupper = round(ARdist[ Int(floor(length(ARdist) * 0.975)) ], digits=4)
            ARlower = round(ARdist[ Int(ceil(length(ARdist) * 0.025)) ], digits=4)
            tpf = round(tpf, digits = 6)
            
            println("Relative risk:\n")
            println("    P($f=$fs)|($g); states=$gs; P = $rrr")
            println("       vs.")
            println("    P($f=$fs)|($g); states=$gso; P = $rro\n")
            printstyled("    Relative Risk: $RR  CI95: $CI95\n\n", color=:cyan)
            println("Absolute risk:\n")
            println("    P($f=$fs)|($g); states=$gs")
            println("       vs.")
            println("    P($f=$fs), all samples; P = $tpf\n")
            printstyled("    Absolute risk: $absrisk  CI95: ($ARlower, $ARupper)\n", color=:cyan)
            println("$('-'^75)")

        end

    end
            
    header = ["Target"]     #create header for multiple states of target
    sck = sum(df_map.feature .== f)
    sc = df_map.state[df_map.feature .== f]

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

    dfp = CSV.read("BN.grain.cpts.out", DataFrame, delim="|", header=header )
    dfp = sort(dfp, [2,3,4], rev=[true,true,true])

    CSV.write("BN_cpts.table", dfp)

    return dfp, df_j, dftf, adjM, mbM, probout  #cpdf, df of used variables, target state freq

end
