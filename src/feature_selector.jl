"""
    dffs, baseprobs  = feature_selector(datatable.csv;  f="target" features=[], ids=[], net="random", netsize=6, relrisk=false, rr_bootstrap=500)

Test every feature as a conditional variable, one-by-one, against a target feature. Iterate over all feature states. Measure the relative change to the baseline probability of the target variable when conditioned on the feature.

The input file is a table of data with obseravation in row and variable/features in columns.  Column 1 must contain the ids and every column labelled. There should be at least 6 variables to test.

All probabilities are network propagated with BNE using the "sm" exact algorithm. Each network will contain the target, the conditional variable, and randomly generated two-state variables. Multi-state random variables can also be selected. Alternatively, variables can be selected randomly from the data set itself (e.g. net="select").

This function returns a data frame (df) of conditional probabilities and baseline probabilities of the target states.  When relrisk=true and net="random", relative and absolute risk ratios with their CI95s are calculated. Relative risk estimates are recorded as the mean value of the bootstrap distribution. It may be necessary to set the netsize to more than the default value if a graph error occurs.

You can use the format_file() function to first create a data frame, a list of ids, and a list of features for input. This data frame can be used for input if the ids are passed to the ids keyword. A subset of features in an input dataframe may be specified in a array (e.g. features=["Age", "Weight"] ) 

Notes: \\

1. Runtimes will increase with netsize. If the frequencies of some variable states are very low, the conditional probability estimates may have high variance. The function implements the SM algorithm and BIC network scoring only.

2. feature_selector() calculates the absolute and relative risk estimates as the mean value of the bootstrap estimates which should be similar to the exact network propagated value from bne if the conditional sample size is sufficient. Feature_selector help to decide what feature to include in a final network. Use bne() to construct the final network and obtain estimates for absolute and relative risks based on that final network.
"""
function feature_selector(filename::Union{String,DataFrame}=""; f::String="",  features::Array=[], ids::Array=[], net::String="random",  netsize::Int64=6, randvarstates::Int64=2, relrisk::Bool=false, rr_bootstrap::Int64=500)

    results=[]
    distributions = Dict()
    
    if typeof(filename) == DataFrame
        if length(unique(filename[!,1])) < size(filename,1)
            if length(ids) == size(filename,1)
                insertcols!(filename, 1, :IDS => ids)
                println("Add ids to the input data frame,")
            elseif length(ids) == 0
                error("Please provide ids for the data frame.")
            else
                error("The ids array must match the data frame.")
            end
        end
        
    elseif isfile(filename)
        dfin = CSV.read(filename, DataFrame, header=true)
        features = names(dfin)[2:end]
        ids = dfin[!,1]
    else
        error("Input must be a complete table with ids and labeled columns or a dataframe. If using a dataframe, you must supply ids and features.")
    end

    if net == "random"
        println("Generating nets with netsize of $netsize")
    elseif   netsize >= (size(filename, 2) - 1) 
        error("\n\nThe total number of test variables must be greater than the requested netsize ($netsize).\n\n")
    elseif netsize < 6
        error("\n\nNetsize should be at least six for random propagation.\n")
    else
    end
    
    if 0 < rr_bootstrap < 100
        error("\n\nUse at least 100 bootstraps to create CI95 estimates, 1000 bootstraps are recommended!\n\n")
    elseif rr_bootstrap >= 100 && relrisk == false
        error("\n\nSet relrisk=true to perform bootstrap analysis.\n\n")
    else
    end
        
    if length(f) > 0
        println("Testing target variable $f conditionally with all provided variables...")
    else
        error("\n\nPlease provide the target variable (e.g. f=\"mytarget\").\n")
    end

    df_all, ids, features_all = format_file(filename) # ids col cropped

    insertcols!(df_all, 1, :IDS => ids)               # replace ids in df
    println("Getting baseline probabilities...")

    fidx = findall(x->x == f, names(df_all))

    if length(fidx) == 0
        error("\n\nDid not find requested target $f in the data!\n")
    else
        fidx = fidx[1]
    end
    
    fcols = [1, fidx]
    fsc = length(unique( df_all[!, fidx ] )) #feature state count

    while length(fcols) <= netsize           #must add variables to make a net!
        rc = rand( 2:size(df_all, 2) )
        if rc == fidx || in(rc, fcols)
        else
            push!(fcols, rc)
        end
    end
    
    fcols = sort(fcols)

    format_file(df_all, datacols=fcols ) #write small BN.data, BN.header 

    cpt, dft, tf, adj, mbM, probout_b = bne("BN.data", "BN.header", f=f);

    probout_b = split(probout_b, "|")
    blist = replace(probout_b[2][2:end-1], " "=>"") #return to Float array
    baseline = parse.(Float64, split( (replace(blist, "\""=>"")), ","))
    
    OUT = open("r.tmp", "w")

    for i in eachindex(features)

        z=[]
        n = features[i]
        givenidx = findall(x->x == n, names(df_all))
        givenidx = givenidx[1]
        gsc = length(unique( df_all[!, givenidx ] )) #given state count

        if fidx == givenidx
            continue
        end
        
        z = [1, fidx, givenidx] 

        println("Processing conditional variable: $n, input column index: $givenidx")

        while length(z) <= netsize          #must add variables to make a net!
            rc = rand( 2:size(df_all, 2) )
            if in(rc, z)
            else
                push!(z, rc)
            end
        end
        
        z = sort(z)
        
        rstates = collect(1:1:randvarstates)

        if net == "random"
            rv = rand( rstates, size(df_all, 1),  (netsize - 2) )
            v = select( df_all, Cols(1, fidx, givenidx) )
            fstates = collect(unique(v[:,2])) 
            v = hcat(v, DataFrame(rv, :auto)) #write sing. var. with random
            format_file(v)                        
        elseif net == "select"        
            format_file(df_all, datacols=z ) #write variable from data         
        else
            error("\n\nValues for net are select or random.\n")
        end
        
        g = [n]
        
        for j in 1:gsc
            gs = []
            push!(gs, string(j))

            ard = []
            rrd = []

            cpt, dft, tf, adj, mbM, probout = bne("BN.data", "BN.header", f=f, g=g, gs=gs);

            if (net == "random") && (netsize < 10) && (netsize > 3) && (relrisk == true)
                
                fstates = sort(fstates)

                for fs in fstates
                    
                    print("Analyzing ==> P(", f, "_", fs, ")|(", join(g,""), "=", join(gs,""), ")\n")

                    cpt, dft, tf, adj, mbM, probout= bne("BN.data", "BN.header"; algo="sm", scoring_method="BIC", f=f, g=g, gs=gs)
                    
                    CI95, PRdist, RR_est = bne_bootstrap("BN.data", "BN.header"; impute=false, algo="sm", scoring_method="BIC", f=f, fs=fs, g=g, gs=gs, bootstrap=0, iterate="", minfreq=0.00, verbose=false, rr_bootstrap=rr_bootstrap, bootstrap_method="resample", boot_plot=true)
                    
                    println("Relative risk ratio: ", RR_est, " ", CI95)
                    push!(rrd, RR_est, CI95[1], CI95[2])

                    fsidx = parse(Int64, fs)
                    baseline_n  = baseline[fsidx]
                    ardist = PRdist ./ baseline_n  #estimates / baseline feature => dist of abs risk

                    AR_est = round(mean(ardist), digits=4)
                    ARlower = round(ardist[ Int(ceil(length(ardist) * 0.025)) ], digits=4)
                    ARupper = round(ardist[ Int(floor(length(ardist) * 0.975)) ], digits=4)

                    push!(ard, AR_est, ARlower, ARupper)
                    println("Absolute risk ratio: $AR_est, ($ARlower, $ARupper)")

                    distributions[i] = PRdist

                end
                
            elseif net == "random" && (netsize < 3 || netsize > 10) && relrisk == true
                error("\n\nDue to graph construction and speed limitation, netsize should be set between 4 and 10 for bootstrapping.\nPlease also make sure there are at least 4 variables in the input data.\n\n")
            end
            
            probout = split(probout, "|")  
            plist = replace(probout[2][2:end-1], " "=>"")
            plist = replace(plist, "\""=>"")
            probout = join([probout[3], probout[4], plist ], ",") 

            if relrisk == true
                rrr = ""
                for i in rrd
                    rrr = rrr * "," * string(i) 
                end
                probout = probout * rrr

                aaa = ""
                for i in ard
                    aaa = aaa * "," * string(i) 
                end

                probout = probout * aaa                

            end

            println(OUT, probout)
            #R"gc()"                       #force R garbage collection

        end

    end

    close(OUT)
    
    tnames = [f * "_" * string(i) for i in 1:fsc]
    tnames2::Vector{String} = []
    
    if relrisk == true
        for i in eachindex(tnames)
            x = tnames[i] * "_RRR"; y = tnames[i] * "_RRR_CI95l"; z = tnames[i] * "_RRR_CI95u";
            xx = tnames[i] * "_ARR"; yy = tnames[i] * "_ARR_CI95l"; zz = tnames[i] * "_ARR_CI95u";
            push!(tnames, x, y, z)
            push!(tnames2, xx, yy, zz)
        end
    end

    tnames = vcat(tnames, tnames2)
    pushfirst!(tnames, "CondVar", "CondVarState")

    df_out = CSV.read("r.tmp", DataFrame, delim=",", header=tnames)
    df_out = sort(df_out, 3)

    for i in eachindex(baseline)
        j = i + 2
        diff_p = round.( (df_out[:,j] .- baseline[i]) .* 100, digits=1)
        dname = tnames[j] * "_" * "diff"
        insertcols!(df_out, dname => diff_p)
    end
    
    rm("r.tmp")
    
    printstyled("Network-propagated baseline probabilities of target states for $f: $baseline\n", color=:green)

    return df_out, baseline, distributions

end
