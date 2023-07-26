"""
    CI95, PRdist = bne_bootstrap("datafile", "headerfile"; impute=false, algo="sm", scoring_method="bic", f="", fs="1",  g=[], gs=[], bootstrap=0, verbose=false, rr_bootstrap=100, bootstrap_method="resample")

    Bootstrap exact Bayesian networks using bnstruct. Inputs are data and header files from makeBNstructfiles.jl. Data must be discrete or quantized if continuous [1,2, ..., N]. 

    This program similar to bne without iteration except that the data set is resampled with replacement or jackknifed (50%) to generate replicates. Relative risk is calculated as the median of the resampled data sets. The CI95 boundaries based on a t-distribution are returned. 

    Note that for each iteration, individuals are randomly resampled, and a new network is remade. Therfore, speed decreases very rapidly as the number of nodes increases.  
"""
function bne_bootstrap(data, header; impute=false, algo=algo, scoring_method=scoring_method, f=f, fs=fs, g=g, gs=gs, bootstrap=0, iterate="", minfreq=0.00, verbose=false, rr_bootstrap=rr_bootstrap, bootstrap_method=boostrap_method, boot_plot=boot_plot, confmeth=confmeth, rrrdenom=rrrdenom)

    if 1 < rr_bootstrap < 100
        printstyled("INFO::For better accuracy, use at least 100 rr_bootstrap replicates.\n", color=:yellow)
    end
    
    if impute == true || bootstrap > 0 || length(iterate) > 0
        error("ERROR::Imputation and iteration with bootstrapping is not yet implemented.")
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
        
        bd = dfcleaned

        # Randomly resample all rows and write to file
        if bootstrap_method == "resample"
           bd = bd[rand(1:nrow(bd), size(bd,1)), :]           
        elseif bootstrap_method == "jackknife"
            nof = Int(floor((size(bd,1)/2)))
            bd = bd[rand(nof, size(bd,1)), :]
        end

        bd = Matrix(bd)

        bd_z = "bd_$z"
        writedlm("bd_$z", bd, " ") #resampled data

        @suppress R"dataset <- BNDataset( $bd_z, $header, starts.from = 1)";   #1-based variables only
        vars = rcopy(R"vars <- variables(dataset)")
        @suppress R"df = read.csv(file = $data, sep = ' '); colnames(df) <- vars" 

        @suppress R"vars <- variables(dataset)";
        @suppress R"net <- learn.network(dataset, algo = $algo, scoring.func = $scoring_method)";
        @suppress net = rcopy(R"net <- learn.network(dataset, algo = $algo, scoring.func = $scoring_method)";)
        vars = rcopy(R"vars <- variables(dataset)")

        @suppress R"dadjm <- net@dag";         
        @suppress R"rownames(dadjm) <- vars; colnames(dadjm) <- vars"; #label DAG matrix
        @suppress R"gn <- as(dadjm, 'graphNEL')";                 # non-boot DAG matrix to graphNEL

        @suppress R"ecpt <- extractCPT(df, gn, smooth = 0.01)";   # extract CPT from data, smooth NaN  
        @suppress R"pt <- compileCPT(ecpt)";
        @suppress R"pnet <- grain(pt, propagate=TRUE)";           # make cpt net, propagated, copy

        if boot_plot == true
            adjM = rcopy(R"dadjm <- net@dag")
            adjM = Int.(adjM)
            plt = plot_network(adjM, headerfile="BN.header", fnode=f, gnodes=g )
            display(plt)
        end
        
        if length(rrrdenom) > 0
            gso = rrrdenom
        else
            gso = replace(gs, "1" => "2", "2" => "1")
        end

        probout, jpout = ConProb(; f=f, fs=fs, g=g, gs=gs, vars=vars, verbose=false, rr_bootstrap=rr_bootstrap);

        probout_o, jpout  = ConProb(; f=f, fs=fs,  g=g, gs=gso, vars=vars, verbose=false, rr_bootstrap=rr_bootstrap);

        RR_all = probout ./ probout_o

        push!(RRdist, RR_all[fs])    #bootstrapped Relative Risk ratios
        push!(PRdist, probout[fs])   #bootstrapped probabilities
        
        @suppress R"rm(list = ls())"    #clear R workspace each run
        rm("bd_$z")

    end

    RRdist = sort(RRdist)
    RRdiststd = std(RRdist)

    if confmeth == "t-dist"
        RR_est = median(RRdist)
        lower, upper = confint(OneSampleTTest(RR_est, RRdiststd, rr_bootstrap), level=0.95)
        if lower < 0
            lower = 0.0
        end        
        CI95 = (round(lower, digits=4), round(upper, digits=4))
    elseif confmeth == "empirical"
        RR_est = median(RRdist)
        upper = RRdist[ Int(floor(length(RRdist) * 0.975)) ]
        lower = RRdist[ Int(ceil(length(RRdist) * 0.025))  ]
        CI95 = (round(lower, digits=4), round(upper, digits=4))
    end

    return CI95, PRdist, RR_est

end
