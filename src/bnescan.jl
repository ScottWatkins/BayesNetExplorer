"""
    ARR, RRR = bnescan(evidence, target)

Calculate P(Y|X), the absolute risk ratio, and the relative risk ratio for a target Y over all evidence variables X = {X1, X2, ... Xn} in an input dataframe. Examime all combinations unless otherwise specified. Return significant absolute and relative risk ratios for any combination. 

        Options:
            ts            target state to evaluate        "all"
            acs           cond. state for all tests       ""
            mincounts     min final, min cond. counts     [0,0]
            outfile       write all results to file       "" bootstraps    number of bootstraps            100
            confmeth      method to get CI95              "t-dist|empirical"
            minup         min pos effect to keep          1.0
            mindown       min neg effect to keep          1.0


"""
function bnescan(e::Union{String,DataFrame}, t::String; ts::String="all", mincounts::Array=[0,0], outfile::String="tmp.out", bootstraps::Int64=100, confmeth="t-dist", acs::String="", minup::Float64=1.0, mindown::Float64=1.0)

    if typeof(e) == String && isfile(e)
        df = CSV.read(e, DataFrame, delim=",", header=true, types=String)
    elseif typeof(e) == DataFrame
        df = string.(e[:,:])
    else
        error("Input must be a csv/tsv file or a dataframe")
    end

    if isfile(outfile)
        rm(outfile)
    end
    
    if minup != 1.0 && mindown != 1.0
        error("\n\nOnly one min filter value allowed per run\n\n")
    end
    
    ac = size(df,2)
    
    if sum(occursin.(t, names(df))) == 1
        println("Analyzing target $t across $ac conditional variables...")
    else
        error("\n\nTarget $t not found in input data.\n\n ")
    end

    
    ats = []
    if ts == "all"
        ats = unique(df[:, Symbol(t)])
        println("Found unique target states: $ts")
    else
        if sum(occursin.(ts, df[:, t] ) ) == 0
            error("\n\nTarget state $ts is not a valid state for $t.\n\n")
        end
        push!(ats, ts)
    end

    for i in names(df)[2:end]

        if i == t
            continue
        end

        if acs == "" 
            cs = unique(df[:, Symbol(i)])  # cond variable states
        else
            cs = [acs]
        end
        
        for j in ats
            for k in cs
                bnemle(df, "P( $t=$j | $i=$k )", mincounts=mincounts, outfile=outfile, bootstraps=bootstraps, confmeth=confmeth)
            end
        end

    end

    labels = string.(split(readlines("header.txt")[1], "\t"))
    dfr = CSV.read(outfile, DataFrame, header=labels);
    dfr = sort(dfr,[2,3])

#println(dfr)
    if minup != 1.0
        dfr_as = dfr[.&(dfr.Estimate .== "ARR (adj)", dfr.Significant .== 1, dfr.DistMean .>= minup), :]
        dfr_rs = dfr[.&(dfr.Estimate .== "RRR (adj)", dfr.Significant .== 1, dfr.DistMean .>= minup), :]
     elseif mindown != 1.0
        dfr_as = dfr[.&(dfr.Estimate .== "ARR (adj)", dfr.Significant .== 1, dfr.DistMean .<= mindown), :]
        dfr_rs = dfr[.&(dfr.Estimate .== "RRR (adj)", dfr.Significant .== 1, dfr.DistMean .<= mindown), :]
    else
        dfr_as = dfr[.&(dfr.Estimate .== "ARR (adj)", dfr.Significant .== 1), :]
        dfr_rs = dfr[.&(dfr.Estimate .== "RRR (adj)", dfr.Significant .== 1), :]
    end
    
    rm("header.txt")

    if outfile == "tmp.out"
        rm(outfile)
        CSV.write("ARR.out", dfr_as)
        CSV.write("RRR.out", dfr_rs)
        println("Key results written to ARR.out and RRR.out\n")
    else
        println("All tests written to $outfile.")
        if isfile("ARR.out")
            rm("ARR.out")
        end
        if isfile("RRR.out")
            rm("RRR.out")
        end
        
    end
    
    
    return dfr_as, dfr_rs
end
