"""
    bnmle(data, "P(A=a|B=b,C=c,...,N=n")

Calculate the maximum likelihood probablity of a target state given any single conditional variable. Input data is csv file or dataframe and a string for the probability query. Only binary conditional variables allow.
"""
function bnmle(data, query; digits=6)

    if occursin(r"P\((.+)\|(.+)\)", query)
        # println(query)
    else
        println("Bad query format -> $query")
    end

    z = queryparser(query)
    # println(z)
    target=""; tstate=""
    for i in eachindex(z[1])
        tt = split(z[1][i], "=")
        target = String.(tt[1])
        tstate = String.(tt[2])
    end

    feature=""; fstate=""
    for i in eachindex(z[2])
        ff = split(z[2][i], "=")
        feature = String.(ff[1])
        fstate = String.(ff[2])
    end
    
    if typeof(data) == DataFrame
        df = data
    elseif isfile(data)
        df = CSV.read(data, DataFrame)
    else
        error("\nCannot find file or dataframe called $data.\n")
    end

    if !in(tstate, df[!, target])
        tv = unique(df[!, target])
        printstyled("\nWARNING: Target state $tstate not found in $target:\nVariables states: $v.\nUsing Laplace correction (k=1) for $target=$tstate", color=:cyan)
    end

    if !in(fstate, df[!, feature])
        fv = unique(df[!, feature])
        error("\nERROR: Feature state $fstate not found in $feature:\nVariables states: $v.\n")
    end
    
    fstate_opp = unique(df[!,feature][findall(df[!, feature] .!= fstate)])

    if length(fstate_opp) == 1
        fstate_opp = fstate_opp[1]
    else
        error("More than two states in selected feature.")
    end
    
    all_c = size(df,1)
    targ_c = sum(df[!, target] .== tstate)
    feat_c = sum(df[!, feature] .== fstate)
    targ_f = targ_c / all_c
    feat_f = feat_c / all_c

    dfc = df[ (df[!, feature] .== fstate), :]
    all_cc = size(dfc, 1)
    targ_cc = sum(dfc[!, target] .== tstate)
    feat_cc = sum(dfc[!, feature] .== fstate)
    targ_fc = targ_cc / all_cc
    feat_fc = feat_cc / all_cc

    dfco = df[ (df[!, feature] .== fstate_opp), :]
    all_cco = size(dfco, 1)
    targ_cco = sum(dfco[!, target] .== tstate)
    feat_cco = sum(dfco[!, feature] .== fstate_opp)
    targ_fco = targ_cco / all_cco
    feat_fco = feat_cco / all_cco

    abs_risk_ratio = targ_fc / targ_f
    rel_risk_ratio = targ_fc / targ_fco

    all_c = size(df,1)
    targ_c = sum(df[!, target] .== tstate)
    feat_c = sum(df[!, feature] .== fstate)
    targ_f = targ_c / all_c
    feat_f = feat_c / all_c
    targ_f = round(targ_f, digits=digits); feat_f = round(feat_f, digits=digits)
    targ_fc = round(targ_fc, digits=digits); feat_fc = round(feat_fc, digits=digits)
    targ_fco = round(targ_fco, digits=digits); feat_fco = round(feat_fco, digits=digits)
    abs_risk_ratio = round(abs_risk_ratio, digits=digits); rel_risk_ratio = round(rel_risk_ratio, digits=digits)

    # With pseudocounts and simple Laplace corrections k=1
    # for zero class problem. That is, given some feature/label
    # the are no target/class instances for at least on feature label state
    
    function Laplace_correction(df, target, tstate, feature, fstate, fstate_opp; k=1)

        states = length(unique(df[!, target]))

        all_c = size(df,1)
        targ_c_lc = sum(df[!, target] .== tstate) + k
        N_lc = all_c + (k * states)
        Ptarg_lc = targ_c_lc / N_lc 
        # println("==>",  targ_c_lc, " ", N_lc, " ", Ptarg_lc)
        # println("\nCorrected probabilty of $target=$tstate = $Ptarg_lc\n")
        
        dfc = df[ (df[!, feature] .== fstate), :]
        all_cc = size(dfc, 1)
        targ_cc = sum(dfc[!, target] .== tstate)
        feat_cc = sum(dfc[!, feature] .== fstate)

        targ_cc = targ_cc + k
        feat_cc = feat_cc + (k * states)  #states of k in target variable 
        Pcond_lc = targ_cc / feat_cc      #matches https://classes.soe.ucsc.edu/cmps140/Winter17/slides/3.pdf
        
        # println("$targ_cc and $feat_cc and $Pcond_lc")
        # println("\nCorrected probabilty of $target=$tstate|$feature=$fstate = $Pcond_lc\n")

        ARR_lc = Pcond_lc / Ptarg_lc
        # println("\nCorrected absolute risk ratio: $ARR_lc\n")

        dfco = df[ (df[!, feature] .== fstate_opp), :]
        # println(size(dfco))
        all_cco = size(dfco, 1)
        targ_cco = sum(dfco[!, target] .== tstate)
        # println(targ_cco)
        feat_cco = sum(dfco[!, feature] .== fstate_opp)
        targ_cco = targ_cco + k
        feat_cco = feat_cco + (k * states)
        Pcondopp_lc = targ_cco / feat_cco
        # println("====>", targ_cco, " ", feat_cco, "    ", Pcondopp_lc)
        RRR_lc = Pcond_lc / Pcondopp_lc

        return (Ptarg_lc, Pcond_lc, ARR_lc, RRR_lc )

    end

    if targ_cc >=  0
        lcvals = Laplace_correction(df, target, tstate, feature, fstate, fstate_opp; k=1)
    end

    lcvals = round.(lcvals, digits=digits)
    # println(lcvals)
    
    line = "-"^70
    println(line)
    println("Base Data Set")
    println(line)
    println(rpad("name=state",30), "\tcount\tfrequency")
    println(rpad("$target=$tstate",30), "\t", targ_c, "\t", targ_f, "  (corrected: ", lcvals[1], ")")
    println(rpad("$feature=$fstate",30), "\t", feat_c, "\t", feat_f)
    println(line)
    println("Conditional Data Set")
    println(line)
    println(rpad("$feature=$fstate", 30), "\t", feat_cc, "\t", feat_fc)
    println(rpad("$feature=$fstate_opp", 30), "\t", feat_cco, "\t", feat_fco)
    println(rpad("$target=$tstate|$feature=$fstate", 30), "\t", targ_cc, "\t", targ_fc, "  (corrected: ", lcvals[2], ")")
    println(line)
    println("Risk Estimates")
    println(line)
    println(rpad("Baseline risk:\nP($target=$tstate)", 57), "$targ_f\n" )
    println(rpad("Conditional:\nP($target=$tstate|$feature=$fstate)", 55), "$targ_fc\n" )
    println(rpad("Abs. risk ratio:\nP($target=$tstate|$feature=$fstate)/P($target=$tstate)", 59), "$abs_risk_ratio", "  (corrected: ", lcvals[3], ")\n" )
    println(rpad("Rel. risk ratio:\nP($target=$tstate|$feature=$fstate)/P($target=$tstate|$feature=$fstate_opp)", 59), "$rel_risk_ratio", "  (corrected: ", lcvals[4], ")\n")

    if size(dfco,1) == 0
        printstyled("WARN: Zero samples found in conditional subset!", color=:magenta )
    elseif 1 <= size(dfc, 1) <= 20
        tiny = size(dfc, 1)
        printstyled("WARNING: Only $tiny samples in conditional subset!\n", color=:magenta )
    end

    if targ_cc == 0
        printstyled("WARNING: Zero target samples in conditional query!\n", color=:magenta )
    elseif 1 <= targ_cc <= 20
        printstyled("WARNING: Only $targ_cc target samples found in conditional subset!\n", color=:yellow )
    end

    println(line)

end
