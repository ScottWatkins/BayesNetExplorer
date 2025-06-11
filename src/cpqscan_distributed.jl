"""
     cpqscan(df, kwargs)

Run cpq on a list of conditional variables and a target outcome.

    Options
        target        the target variable                       ""
        targetstate   the target state                          ""
        clist         vec array of conditional variables        []
        cstates       vec array of conditional states           []
                      (defaults to "1" unless specified)
        outfile       outfile name                              scanout.csv
        mincounts     final and conditional minimum counts      [0,0]
                      (filtered at end)
        decode        add conditional state terms from file     ""
                      ("/path/file.tsv,cstatecol,altnamecol")
        cores         number of cores                           4

"""
function cpqscan(df::DataFrame; target::String="", targetstate::String="", clist::Array=[], cstates::Array=[], outfile::String="scanout.csv", mincounts::Array=[0,0], decode::String="", cores=4)

    addprocs(cores)
    @everywhere 
    
    if isfile(outfile)
        rm(outfile)
    end
    
    
    t = target
    ts = length(targetstate) == 0 ? "1"  : targetstate    
    println("Target variable is $t ...\nTarget state set to $ts ...")

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
        println(c, "  ", cs)
        @suppress x = cpq(df, "P($t=$ts|$c=$cs)", bootstraps=100, outfile=outfile, mincounts=[0,0], binomial=true);
    end
    println("\ndone.")

    h = string.(split(read("header.txt", String), "\t"))
    h[end] = chomp(h[end])

    dfa = CSV.read(outfile, DataFrame, header=h)

    rm("scanout.csv")
    rm("header.txt")

    if length(decode) > 0
        f, ncol, nncol = string.(split(decode, ","))

        ncol = parse(Int, ncol)
        nncol = parse(Int, nncol)
        
        lup = hashfile(f, ncol, [nncol]) #hashfile modified to match colnames

        rc = []

        for i in clist
            r = get(lup, i, i)
            push!(rc, r)
        end

        insertcols!(dfa, 2, "CondTerms" => rc )
        dfa = sort(dfa[:, 1:11], 7, rev=true)

    else
        
        dfa = sort(dfa[:, 1:10], 7, rev=true)

    end

    if mincounts != [0,0]
        dfa = dfa[.&(dfa.Count .> mincounts[1], dfa.CondCount .> mincounts[2]), :]
    end
    
    CSV.write(outfile, dfa)
    println("Results written to $outfile.")

    return dfa

end
