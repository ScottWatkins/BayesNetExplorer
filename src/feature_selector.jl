"""
    df = feature_selector(filename;  f="target", features=[], ids=[], netsize=7)

Test each feature (variable) in a data set, one-by-one, against the target variable. Measure the relative change to the baseline probability of the target variable when conditioned on the feature. Probabilities are network propagated with BNE.

This function returns a data frame (df) of conditional probabilities and the change relative to the target variable baseline probability over all target and conditional states.  

Notes: please use the format_file() function to create the BN.data file, the BN.header files, and list of features for input. If a crash occurs due to graphNel, try increasing the netsize. Runtimes will increase with increasing netsize.

"""
function feature_selector(filename::Union{String,DataFrame}=""; f::String="",  features::Array=[], ids::Array=[], netsize::Int64=7)

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
    end
    
    results=[]
    
    if length(f) > 0
        println("Testing target variable $f conditionally with all provided variables...")
    else
        error("\n\nPlease provide the target variable (e.g. f=\"mytarget\").\n")
    end

    df_all, ids, features_all = format_file(filename) # ids col cropped

    insertcols!(df_all, 1, :IDS => ids)               # replace ids in df
    println("Getting baseline probabilities...")

    fidx = findall(x->x == f, names(df_all))
    fidx = fidx[1]
    fcols = [1, fidx]
    fsc = length(unique( df_all[!, fidx ] )) #feature state count

    while length(fcols) <= netsize          #must add variables to make a net!
        rc = rand( 2:size(df_all, 2) )
        if rc == fidx || in(rc, fcols)
        else
            push!(fcols, rc)
        end
    end
    
    fcols = sort(fcols)

    format_file(df_all, datacols=fcols ) #rewrite small input BN.data, BN.header 

    cpt, dft, tf, adj, mbM, probout_b = bne("BN.data", "BN.header", f=f);
    probout_b = split(probout_b, "|")
    blist = replace(probout_b[2][2:end-1], " "=>"") #ugh, return to Float array
    baseline = parse.(Float64, split( (replace(blist, "\""=>"")), ","))
    
    OUT = open("r.tmp", "w")

    for i in eachindex(features)

        z=[]
        n = features[i]
        givenidx = findall(x->x == n, names(df_all))
        givenidx = givenidx[1]
        gsc = length(unique( df_all[!, givenidx ] )) #feature state count

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

        format_file(df_all, datacols=z ) #rewrite a single variable BN.data, BN.header
          
        g = [n]

        for j in 1:gsc
            gs = []
            push!(gs, string(j))

            cpt, dft, tf, adj, mbM, probout = bne("BN.data", "BN.header", f=f, g=g, gs=gs);

            probout = split(probout, "|")  
            plist = replace(probout[2][2:end-1], " "=>"") #weird! WTF
            plist = replace(plist, "\""=>"")
            probout = join([probout[3], probout[4], plist ], ",") 
            
            println(OUT, probout)
            #R"gc()"                       #force R garbage collection
        end

    end

    close(OUT)

    
    tnames = [f * "_" * string(i) for i in 1:fsc]
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

    printstyled("Baseline probabilities for $f target states: $baseline\n", color=:green)

    return df_out

end
