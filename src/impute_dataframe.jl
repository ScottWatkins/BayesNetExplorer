"""
    dfi = impute_dataframe(df; boolout=false)

Imputes missing values in a dataframe. Input should be Float64 or Int64.
Input must have fully labeled observations and variables!

Options:
    boolout    output will be true/false instead of 0/1
    k          number of nearest neighbors [10]
"""
function impute_dataframe(df; boolout=false, k=10)

    dfnn = df[!, 2:end]
    dfmf = Matrix{Union{Float64,Missing}}(dfnn)          #Int to Float64
    dfmi = Impute.knn(dfmf; dims=:cols, k=k)             #k-nearest neighbors
    dfmi = DataFrame(Matrix{Int}(floor.(dfmi)), :auto)   #back to Int
    rename!(dfmi,  names(df)[2:end])                     #return col names
    insertcols!(dfmi, 1, :BlindedId => df[:,1])          #return rownames

    if boolout == true
        dfmi[!,2:end] = Bool.(dfmi[!,2:end])
    end

    CSV.write("BN.imputed.csv", dfmi, delim=",")          #save to disk
    println("Wrote BN.imputed.csv to disk.")

    return dfmi

end
