"""
    correlation_analysis(df_r, num_samples; kwargs)

Analyze a square matrix of correlations in an input dataframe (row and columns labeled). Calculate p-values and flag the colinear pairs. Make a heat map showing correlation and pvalues.

    Options
        pmax    max p-values        [0.05]
        rmin    min cor for flag    [0.7]
        title   plot title          ""
        cbtitle colorbar title      "Correlation"

    Input Example:
         Row │ Features  A         B          C          D         E        
             │ Symbol    Float64   Float64    Float64    Float64   Float64  
        ─────┼──────────────────────────────────────────────────────────────
           1 │ A         1.0       0.270071   0.376995   0.822757  0.287241
           2 │ B         0.270071  1.0        0.0229777  0.358427  0.571194
           3 │ C         0.376995  0.0229777  1.0        0.501568  0.985098
           4 │ D         0.822757  0.358427   0.501568   1.0       0.21191
           5 │ E         0.287241  0.571194   0.985098   0.21191   1.0

"""
function correlation_analyzer(df_r::DataFrame, num_samples::Int64; pmax::Float64=0.05, rmin::Float64=0.7, tails::Int64=2, title::String="", cbtitle::String="Correlation", fontsize=6)

    
    cbtitle = "\n" * cbtitle
    n = names(df_r)
    M = Matrix(df_r[:, 2:end])
    size(M,1) == size(M,2) ? println("Input contains symmetric correlation matrix...") : error("Correlation matrix is not symmetric!")
    
    print("Calculating p-values...")
    M_pvals = cor_p.(abs.(M), num_samples, tails)
    println("done.")

    df_p = DataFrame(M_pvals, :auto)
    insertcols!(df_p, 1,  :Features => names(df_r)[2:end] )
    fnames = Dict(zip(["x" .* string(i) for i in 1:(size(df_r,2))], names(df_r)[2:end]))
    rename!(df_p, fnames)
    vnames = names(df_r)[2:end]

    nada = fill(NaN, size(M))
    
    function merge_tril_triu(L, U)
        MM = ones(size(L))
        for i in 1:size(MM,1)
            for j in 1:size(MM,2)
                if i>j
                    MM[i,j] = L[i,j]
                elseif i<j
                    MM[i,j] = U[i,j]
                end
            end
        end
        return MM
    end
    
    MM = merge_tril_triu(nada, M) #set pvals to NaN to prevent color map
    
    plt_r = heatmap(1:size(MM,1), 1:size(MM,2), MM,
                    yflip=true,
                    xticks=(1:size(MM,1), vnames),
                    xrotation=45,
                    yticks=(1:size(MM,2), vnames),
                    title=title,
                    c = cgrad(:heat),
                    tickfontsize=10,
                    tickorientation=:out,
                    colorbar_title=cbtitle,
                    right_margin = 5Plots.mm,
                    #top_margin = 5Plots.mm,
                    bottom_margin = 3Plots.mm,
                    grid=:off,
                    framestyle=:box
                  
                    )

    glines = collect(1:1:(size(MM,1) - 1)) .+ 0.50
    vline!(glines, label="", color=:grey50, alpha=0.5)
    hline!(glines, label="", color=:grey50, alpha=0.5)
    # Note j and i are flipped to match flip_y above!!
    MM = merge_tril_triu(M_pvals, M)    #remake pval, correlation matrix

    # Use sprintf to format to decimal and convert to string values; no rounding after this
    
    M2 = [ @sprintf("%-3.3f", MM[i,j]) for j in 1:size(MM,2) for i in 1:size(MM,1) ]
    MM = reshape(M2, size(MM,1), size(MM,2))

    nrow, ncol = size(MM)

    ann = [(j,i, text((MM[i,j]), fontsize, :black, :bold, :center)) for i in 1:nrow for j in 1:ncol]

    annotate!(ann, linecolor=:blue, linewidth=2)

    display(plt_r)

    println(df_r)
    println(df_p)
    
    OUT = open("CorrelatedPairs.csv", "w")
    println(OUT, "Pair,Correlation,P-value,Colinear")

    for i in 1:size(M,1)
        for j in 1:size(M,2)
            if i == j
                continue
            elseif M_pvals[i,j] <= pmax && M[i,j] <= rmin
                println(OUT, vnames[i], "-", vnames[j], ",", M[i,j], ",", M_pvals[i,j], ",No")
            elseif M_pvals[i,j] <= pmax && M[i,j] > rmin
                println(OUT, vnames[i], "-", vnames[j], ",", M[i,j], ",", M_pvals[i,j], ",YES")
            else
                print(".")
                #println(vnames[i], " and ",  vnames[j], " are not significant.")
            end
        end
    end
    close(OUT)
    println("\nCorrelated pairs written to file: CorrelatedPairs.csv")
    return df_p, plt_r 

end
