"""
    getcounts(df_cpt, dfc, target_state)

Input is conditional probabilities from bne, translation_table from bne, and the RRcalculator input of target_state (e.g. wetgrass_true) where wetgrass is the target and "true" is the state.  

Use this function in RRcalculator to get conditional counts and final counts.
"""
function getcounts(df_cpt, dfc, target_state)

    t = split(target_state, "_")
    fs = string(t[end])

    if occursin(r"(2)|(true)|(yes)"i, fs)
        fs = "2"
    elseif occursin(r"(1)(false)|(no)"i, fs)
        fs = "1"
    else
        error("If X is a target, the target-state for true can be X_true, X_yes, or X_2 (case insensitive).\n Only two-state variables are allowed for risk calculations")
    end
    
    #println("^^", fs)
    
    f = join(t[1:end-1], "_")

    tcounts = Int64[]
    condcounts = Int64[]

    for i in 1:size(df_cpt,1)

        cv = split(df_cpt[i,7], ",")
        cs = split(df_cpt[i,8], ",")
        
        q = ""

        for i in 1:length(cv)
	    q = q * "dfc." * cv[i] * ".==" * cs[i] * ","
        end
            
        q = "dfc[ .&(" * q[1:end-1] * "), :]"

        dcall = Meta.parse(q)
        global dfsubset = eval(dcall)

        p = "dfsubset[dfsubset." * f *  " .== " * fs * " , :]"
        tcall = Meta.parse(p)
#        println("===>",tcall)
        fcount = eval(tcall)
        
        push!(condcounts,  size(dfsubset,1))
        push!(tcounts, size(fcount,1))
    end

    return condcounts, tcounts

end
