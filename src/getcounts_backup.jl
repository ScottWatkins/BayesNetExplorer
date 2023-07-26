function getcounts(df_cpt, dfc)

    rawcounts = Int64[]

    for i in 1:size(df_cpt,1)

        cv = split(df_cpt[i,7], ",")
        cs = split(df_cpt[i,8], ",")
        
        q = ""

        for i in 1:length(cv)
	    q = q * "dfc." * cv[i] * ".==" * cs[i] * ","
        end
            
        q = "dfc[ .&(" * q[1:end-1] * "), :]"

        dcall = Meta.parse(q)
        dfsubset = eval(dcall)
        push!(rawcounts,  size(dfsubset,1))
        
    end

    return rawcounts

end
