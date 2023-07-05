"""
     pairwise_correlations(df::DataFrame; method="Pearson")

Calculate the pairwise correlations among variable in a dataframe.

    Options
        method   algorithm to use    Pearson|Spearman|Kendall

"""
function pairwise_correlations(df::DataFrame; method="Pearson")

    n = names(df)
    M = Matrix(df)

    if method == "Pearson"
        C = cor(M)
    elseif method == "Spearman"
        C = corspearman(M)
    elseif method == "Kendall"
        C = corkendall(M)
    else
        println("Please use Pearson, Spearman, or Kendall for method.")
    end

    pwc = DataFrame(C, :auto)
    vnames = Dict(zip([Symbol.("x" .* string(i)) for i in 1:(size(M,2))], n))

    pwc =  rename(pwc, vnames)
    insertcols!(pwc, 1, :Features => n)
    return pwc

end
