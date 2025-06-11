"""
     cpqscan(df, kwargs)

Run cpq using list of conditional variables and a target outcome.

_Fullscan_
1. User simply specifies a single target variable and its targetstate. All other variables are tested for pairwise conditional association with the target.

cpqscan(dft, target="Smoking", targetstate="1")

_Custom scan_
1. Make an array of selected variables to test \\
   e.g., cols = names(df) \\
   selectedvars = cols[occursin.("ICD10", cols)] \\

cpqscan(dft, target="Smoking", targetstate="1", clist=selectedvars, decode="/path/terms_names.tsv,1,4")  

Hint: for large data sets, set bootstraps=1 to increase speed then retest select variables with more bootstraps.

    Options
        fullscan      scan all against target, targetstate      "true"

        target        the target variable                       ""
        targetstate   the target state                          ""

        clist         vec array of conditional variables        []
        cstates       vec array of conditional states           []
                      (defaults to "1" unless specified)
        outfile       outfile name                              scanout.csv
        mincounts     final and conditional minimum counts      [0,0]

        decode        add conditional state terms from file     ""
                      ("/path/file.tsv,codecolumn,namecolumn")

"""
function cpqscan(df::DataFrame; target::String="", targetstate::String="", clist::Array=[], cstates::Array=[], outfile::String="scanout.csv", mincounts::Array=[0,0], decode::String="", bootstraps::Int64=1, fullscan::Bool=true)

    if isfile(outfile)
        rm(outfile)
    end

    scor = []

    t = target
    ts = targetstate
    println("Target variable is $t set to targetstate $ts ...")
    cnames = names(df)

    tidx = findall(occursin.(t, cnames))
    
    if length(tidx) == 1     # check target and get index
        println("Found requested target $t at column $tidx ...")
    else
        error("\n\nRequested target $t not found in dataframe.\n\n")
    end

    
    if fullscan == true && length(clist) == 0

        println("Performing a full scan ...")
        
        for i in 2:size(df,2)

            i % 10 == 0 ? print("...", i) : nothing

            if i == tidx[1]

                println("Skipping column $i ...")
                continue

            else

                dfs = df[!, [1, tidx[1], i]]            # extract pairs
                sc = corspearman(parse.(Int, dfs[!,2]), parse.(Int, dfs[!,3]))
                push!(scor, sc)
                c = cnames[i]
                cs = "1"

                @suppress x = cpq(df, "P($t=$ts|$c=$cs)", bootstraps=bootstraps, outfile=outfile, mincounts=[0,0], binomial=true, digits=6);
                           
            end
            
        end

        println(" done.")
    
    elseif fullscan == false && length(clist) > 0

#        t = target
#        ts = length(targetstate) == 0 ? "1" : targetstate    
#        println("Target variable is $t ...\nTarget state set to $ts ...")

        lc = length(clist)
        cstates = length(cstates) == 0 ? repeat(["1"], lc) : cstates
        lcs = length(cstates)
        println("Conditional variables: $lc ...\nConditional states: $lcs ...")
        
        if lc != length(cstates)
            error("\n\nNumber of conditional variables must match number of conditional states.\n\n")
        end
        
        print("Performing scan...")
        for i in eachindex(clist)
            i % 100 == 0 ? print("..", i) : nothing
            c = clist[i]
            cs = cstates[i]
            typeof(c) == Symbol ? c = string(c) : nothing
            cidx = findall(occursin.(c, names(df))) # one at a time
            print(".")

            dfs = df[!, [1, tidx[1], cidx[1]]]      # extract target-conditional pairs
            sc = corspearman(parse.(Int, dfs[!,2]), parse.(Int, dfs[!,3]))
            push!(scor, sc)

            @suppress x = cpq(df, "P($t=$ts|$c=$cs)", bootstraps=bootstraps, outfile=outfile, mincounts=[0,0], binomial=true, digits=6);

        end        
        println("\ndone.")

    elseif (fullscan == true && length(clist) > 0) || (fullscan == false && length(clist) == 0)
        error("\n\nPlease do a fullscan (default)\nOR provide a vector of variables (clist=[:var1, :var2]) and set fullscan=false.\n\n")
    else
        error("\n\nAn unspecified error occurred.\n\n")
    end
    
    h = string.(split(read("header.txt", String), "\t"))
    h[end] = chomp(h[end])

    dfa = CSV.read(outfile, DataFrame, header=h)
    
    rm("header.txt")

    if length(decode) > 0
        
        decode = replace(decode, " " => "")

        if count(r",|\t", decode) == 2
            println("Reading decode file...")
        else
            error("\n\nThe comma (or tab) delimited decode format should be\n\"filename,codecolumn,termcolumn\"\n(e.g, decode=\"terms.csv,1,4\")\n\n")
        end
        
        f, ncol, nncol = string.(split(decode, r",|\t"))
            
        ncol = parse(Int, ncol)
        nncol = parse(Int, nncol)
        
        lup = hashfile(f, ncol, [nncol])  #hashfile modified to normalize keys
        rc = []

        if length(clist) > 0
            cinfo = clist
        else
            cinfo = names(df)[2:end]
            println("How many additional variable did you create and add to the dataframe?")
            e = parse(Int, readline()) -1
            for i in 1:e
                popfirst!(cinfo)
            end
        end

        for i in cinfo
            r = get(lup, i, "NA")
            push!(rc, r)
        end

        insertcols!(dfa, "Terms" => rc )
        
    else
    end

    insertcols!(dfa, :Correlation => scor)    
    impARR = dfa.cpq_ARR_est .* (dfa.CondCount ./ size(df,1))
    insertcols!(dfa, :ARRimpact => impARR)
    
    if mincounts != [0,0]
        dfa = dfa[.&(dfa.Count .> mincounts[1], dfa.CondCount .> mincounts[2]), :]
    end

    dfa = sort(dfa, 11, rev=false)
    
    dfa = select(dfa, Not(11:16))  # remove cols with missings
    dfa = sort(dfa, [10,3,2], rev=[false,true,true])
    
    CSV.write(outfile, dfa)
    println("Results written to $outfile.")
    println("Important: run at least 100 bootstraps to get final ARR/RRR estimates!")

    return dfa
        
end
