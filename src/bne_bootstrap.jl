"""
    CI95, PRdist = bne_bootstrap("datafile", "headerfile"; impute=false, algo="sm", scoring_method="bic", f="", fs="1",  g=[], gs=[], bootstrap=0, verbose=false, rr_bootstrap=100, bootstrap_method="resample")

    Bootstrap exact Bayesian networks using bnstruct. Inputs are data and header files from makeBNstructfiles.jl. Data must be discrete or quantized if continuous [1,2, ..., N]. 

    This program similar to bne without iteration except that the data set is either resampled with replacement or jackknifed (50%) to generate replicates. Relative risk is calculated for each resampled data set and the CI95 boundaries are returned. The distribution of target probabilities are also returned.

    Network learning methods:\\
        sm (Silander-Myllymaki exact search)\\
        hc (Hill Climbing method)\\
        mmhc (Max-Min Hill-Climbing heuristic, default)\\
        sem (Structural Expectation-Maximization, uses EM for imputation, can set EM params)\\

    Scoring methods:\\
        BDeu (Bayesian-Dirichlet equivalent uniform, default)\\
        AIC (Akaike Information Criterion)\\
        BIC (Bayesian Information Criterion)\\

    Other options:
        minfreq      min variable state frequency cutoff [0.00]
        impute       impute with bnstruct [false]
    
    Probability options with input example:
        f            Target feature (required)           "F1"
        g            Given/Evidence features            ["F2", "F3"]
        gs           Given/Evidence  states             ["2",   "1"]

    Examples:

    1. Simple query for a target and two evidence features and their states:
    bne("BN.data", "BN.header", algo="sm", scoring_method="BIC", f="A", g=["B", "C", ], gs=["2","1"] )

    2. Query feature for 3 specified nodes. All states for all nodes examined.
    bne("BN.data", "BN.header", algo="hc", scoring_method="AIC", f="A", g=["B", "C", "D" ], iterate="nodes" )

    2. Query feature for all nodes. All states for all nodes examined.
    bne("BN.data", "BN.header", algo="sm", scoring_method="BDeu", f="A", iterate="all" )


"""
function bne_bootstrap(data, header; impute=false, algo=algo, scoring_method=scoring_method, f=f, fs=fs, g=[], gs=[], bootstrap=0, iterate="", minfreq=0.00, verbose=false, rr_bootstrap=rr_bootstrap, bootstrap_method=boostrap_method)

    if rr_bootstrap < 10
        error("For accuracy, please use at least 10 rr_bootstrap replicates.")
    end
    
    if impute == true || bootstrap > 0 || length(iterate) > 0
        error("Imputation, iteration, and network bootstrapping not yet implemented.")
    end

    R"library(bnstruct)";  R"library(qgraph)"; R"library(gRain)"
    
    vardat = readlines("$header")
    headdat = replace(readline("$header"), "\"" => "")
    headdat = String.(split(headdat, r"\s+"))
    numstates = parse.(Int64, split(vardat[2]))
    types = split(vardat[3])
    
        dfcleaned, df_freq = clean_dfvars_by_frequency("$data", minfreq; nolabels=true, delim=" ", header=headdat)

    
    RRdist = Float64[]   
    PRdist = Float64[]
    fs = parse(Int, fs)

    print("Bootstrapping: ")        

    for z in 1:rr_bootstrap
        
        if z % 10 == 0
            print(z)
        end

        z == rr_bootstrap ? println() : print(".")

        bd = Matrix(dfcleaned)
        
        if bootstrap_method == "resample"
            bd = bd[sample(axes(bd, 1), size(bd,1); replace = true, ordered = false), :]
        elseif bootstrap_method == "jackknife"
            nof = Int(size(bd,1)/2)
            bd = bd[sample(axes(bd, 1), nof), :]
        end

        writedlm("bd_$z", bd, " ")

        bd_z = "bd_$z"
        @suppress R"dataset <- BNDataset( $bd_z, $header, starts.from = 1)";  #must use 1-based variables
        vars = rcopy(R"vars <- variables(dataset)");

        if impute == true
            @suppress R"dataset <- impute(dataset)";
        end

        if impute == true && bootstrap > 0
            @suppress R"dataset <- bootstrap(dataset, num.boots = $bootstrap, imputation = $impute)";
            @suppress R"net <- learn.network(dataset, algo = $algo, scoring.func = $scoring_method, use.imputed.data = $impute)";
        elseif impute == false && bootstrap > 0
            @suppress R"dataset <- bootstrap(dataset, num.boots = $bootstrap)";
            @suppress R"net <- learn.network(dataset, algo = $algo, scoring.func = $scoring_method)";
        elseif impute == true && bootstrap == 0
            @suppress R"net <- learn.network(dataset, algo = $algo, scoring.func = $scoring_method, use.imputed.data = $impute)";
        elseif impute == false && bootstrap == 0
            @suppress R"net <- learn.network(dataset, algo = $algo, scoring.func = $scoring_method)";
        end

        boots = rcopy(R"b <- has.boots(dataset)")
        bootcount = rcopy(Int.(R"bc <- num.boots(dataset)"))
        
        if boots == true
            @suppress R"wpnet  <- wpdag.from.dag(net)";  #wp net object
            @suppress R"wpadjm <- wpnet@wpdag";          #bootstrapped adj. matrix
            @suppress R"dadjm  <- wpnet@dag";            #non-bootstraped adj. matrix
            @suppress R"rownames(wpadjm) <- vars; colnames(wpadjm) <- vars"; #label boot matrix
            @suppress R"rownames(dadjm) <- vars; colnames(dadjm) <- vars";   #label non-boot matrix
        elseif boots == false
            @suppress R"dadjm <- net@dag";                #non-bootstrapped DAG adj matrix
            @suppress R"rownames(dadjm) <- vars; colnames(dadjm) <- vars"; #label DAG matrix
        end

        @suppress R"df = read.csv(file = $data, sep = ' ')";      # get data table
        @suppress R"colnames(df) <- vars";                        # label data table
        @suppress R"gn <- as(dadjm, 'graphNEL')";                 # non-boot DAG matrix to graphNEL
        @suppress R"ecpt <- extractCPT(df, gn, smooth = 0.01)";   # extract CPT from data, smooth NaN  
        @suppress R"pt <- compileCPT(ecpt)";
        @suppress R"pnet <- grain(pt, propagate=T)";              # make cpt net, propagated, copy

        gso = replace(gs, "1" => "2", "2" => "1")

        check = unique(gso)
        if length(check) > 2
            error("Only 2-state variables are currently allowed.")
        end
        
        probout = ConProb(; f=f, g=g, gs=gs, vars=vars, verbose=false, rr_bootstrap=rr_bootstrap);

        probout_o = ConProb(; f=f, g=g, gs=gso, vars=vars, verbose=false, rr_bootstrap=rr_bootstrap);
        
        RR_all = probout ./ probout_o

        push!(RRdist, RR_all[fs])
        push!(PRdist, probout[fs])
        
        R"rm(list = ls())"    #clear R workspace each run

    end
    
    RRdist = sort(RRdist)
    upper = RRdist[ Int(floor(length(RRdist) * 0.975)) ]
    lower = RRdist[ Int(ceil(length(RRdist) * 0.025)) ]
    
    CI95 = (round(lower, digits=4), round(upper, digits=4))
    RR_est = round(mean(RRdist), digits=4)
    run(`bash -c "rm -f bd_*"`)
    
    return CI95, PRdist, RR_est

end

