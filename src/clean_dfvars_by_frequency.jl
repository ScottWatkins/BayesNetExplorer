"""

    dfcleaned, dfreq = clean_dfvars_by_frequency(df::Union{DataFrame,String}, minfreqpercent::Float64, nolabels=false, header=true)

Remove variables from a dataframe if a variable as a state with a frequency less than minfreq.

"""    
function clean_dfvars_by_frequency(data::Union{DataFrame,String}, minfreq::Float64; nolabels=false, header=true, delim=",")

    minfreq = minfreq*100; labels=[]

    if typeof(data) == String && isfile(data)

        dfcheck = CSV.read(data, DataFrame, delim=delim, header=header)

        if nolabels == true
            dfcheck = dfcheck
        else
            labels = dfcheck[:,1]
            dfcheck = dfcheck[:, 2:end]
        end        

    else

        if nolabels == true
            dfcheck = data
        else
            labels = data[:,1]
            dfcheck = data[:, 2:end]
        end
    end
    
    rmcols = []; df_freq = DataFrame();

    for i in 1:size(dfcheck,2)

        tvar = names(dfcheck)[i]
        ftab = transform!(combine(groupby(dfcheck, [i]), nrow => :count), (:count => (x -> (x/sum(x))*100) => :percent))
        ftab = insertcols!(ftab, 1, :feature => [ names(ftab)[1] for i in 1:size(ftab,1)] )
        rename!(ftab, Symbol(names(ftab)[2]) => :numstate)

        if (sum(ftab.percent .< minfreq)) > 0
            push!(rmcols, Symbol(tvar) )
            printstyled("Warning: $tvar has a variable state below the requested minimum frequency ($minfreq%)!:\n", ftab, "\n", color=:yellow)
        else
            df_freq = vcat(df_freq, ftab)
        end

    end

    dfcleaned = select(dfcheck, Not(rmcols))

    return dfcleaned, df_freq, labels

end
