"""
    pc = bne_cor(table; method="pearson")

Read the formatted dataframe of observations and variables and return all pairwise correlations in a dataframe.

Input: Observations in rows, variables in columns. The dataframe is the bne input dataframe from the format_data() function and does not contain sample names

Methods: pearson, spearman, kendall
"""
function bne_cor(table; method="pearson")

    if typeof(table) == String
        df = CSV.read(table, DataFrame, comment="#")
    elseif typeof(table) == DataFrame
        df = table
    else
        error("\nInput must be a fully label file or dataframe.")
    end
    
    M = Matrix(df[:,1:end]) 
    
    if typeof(M[2,2]) <: AbstractString
        println("Detected string data. Converting to numeric data...")
        if sum(occursin.(r"(yes)||(no)"i, M[:,1])) > 0
            M = parse.(Int64, replace.(M, r"no"i => "0", r"yes"i =>"1"))
        else
            error("Only binary yes/no string data can be automatically converted")
        end
    else
        println("Getting pairwise correlations...")
    end

    println("Input data sample... \n", M[1:3, 1:3])
    
    if method == "pearson"
        cM = cor(M)
    elseif method == "spearman"
        cM = corspearman(M)
    elseif method == "kendall"
        cM = corkendall(M)
    else
        error("\nMethods are pearson, spearman, or kendall\n\n")
    end

    eM = triu(cM, 1)
    s = round.(sort(vcat(eM...)), digits=4)
    println("Bottom and top 3 values ...")
    println( join(s[1:3], ","), " ... ", join(s[end-3:end], ",") )
    println("Range for correlations: ", round.( extrema(eM) , digits=4) )
    
    cn = names(df)[1:end]
    cM = DataFrame(cM, cn)
    insertcols!(cM, 1, :Features => cn)

    return cM
    
end
