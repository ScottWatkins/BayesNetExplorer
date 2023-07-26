"""
    bnenbc(data, "P(Y=y1|X1=x1, X2=x2)")

A naive Bayes classifier that provides estimates of the posterior probability for an outcome (Y) given some body of evidence (X) that contains categorical variables, X = {x1, x2, ..., xn}. 

Inputs: a dataframe with ids in column 1 and the probability query.

If a second query is provided, the relative risk of the posterior probability of query1 vs. query2 is calculated. Class conditional probabilities with zero observations are Laplace corrected. Variables with >2 states are handled by a 2-state, one vs. the rest approach.

    Keyword options:
        digits           Round to n digits                      [4]    
        outfile          Output file                            ""  
        bootstraps       Number of replicates                   [100]  
        k                Laplace correction factor              [1]  
                         (Set k=2 for Wilson midpoint)
        q2               Custom rel. risk query                 "query"
                         [typically, P(Y|X1=No,X2=No...)]
"""
function bnenbc(data::Union{String,DataFrame}, query::String; q2::String="", digits::Int64=6, k::Number=1, outfile::String="", bootstraps::Int64=0, mincounts::Array=[1,20], outfile2::String="" )

 
    if typeof(data) == DataFrame
        df = data
    else
        error("Input data variable must be dataframe.")
    end
    
    digits > 8 ? error("Significant digits maximum is 8 digits.") : "";

    if bootstraps > 0
        if length(q2) < 6
            error("\n\nA second query, q2=P(Y...|X...), must be specified\nfor bootstrapping relative risk.\n\n")
        end
    end


    line = "-"^78
    println(line)
    
    function getinput(input_query)
        
        z = queryparser(input_query)
        println(line)
    
        ta = split(z[1][1], "=")
        t  = string(ta[1])
        ts = string(ta[2])       # t = target name, ts = target state
    
        f = Dict()                 # Dict of all conditional features and states

        for i in z[2]
            x =  split(i, "=")
            f[x[1]] = x[2] 
        end
        
        return t, ts, f

    end

    q1t, q1ts, q1f = getinput(query)    

    function getcounts(dfin::DataFrame, target::String, targetstate::String, features::Dict, k::Int64)
        df = dfin; t = target; ts = targetstate; f=features; k=k;
        rc   = Dict()    # raw counts
        cf   = Dict()    # corrected freq
        rc_o = Dict()    # raw counts !Y
        cf_o = Dict()    # corrected freq !Y
        
        in(ts, df[!, t]) ? "" : error("\n\n$t=$ts not found in data.\n\n")

        rows = size(df,1)
        tc = sum(df[!,t] .== ts) # target counts given target state
        tf = tc / rows           # target state freq

        tc_o = sum(df[!,t] .!= ts)   # target counts opp
        tf_o = tc_o / rows           # target state freq opp

        println()
        println(line)
        println(rpad("Query", 45), rpad("Counts", 15),  "Frequency")
        println(line)
        println("Prior probabilities [P(Y)]\n")
        println(rpad("    $t=$ts", 45), rpad("$tc / $rows", 15), round(tf, digits=digits))
        println(rpad("    $t=!$ts", 45), rpad("$tc_o / $rows", 15), round(tf_o, digits=digits), "\n")
        println("Class conditional probabilities [P(X|Y)]\n")
        
        for i in keys(f)         # class cond. probs

            fs = f[i]   
            in(fs, df[!, i]) ? "" : error("\n\n$i=$fs not found in data.\n\n")
            c = sum(df[!, i] .== fs .&& df[!, t] .== ts)   # counts of Xi|Y
            d = sum(df[!, t] .== ts)                       # counts of Y where t=ts

            c_o = sum(df[!, i] .== fs .&& df[!, t] .!= ts)   # counts of Xi|Y
            d_o = sum(df[!, t] .!= ts)                       # counts of Y where t=ts
            
            if c == 0                       
                n = length(unique(df[!, i]))
                c = k + c            # class cond. prob correction for numerator
                d = (k * n) + d      # and denominator
            end
            
            if c_o == 0
                n = length(unique(df[!, i]))
                c_o = k + c_o            # class cond. prob correction for numerator
                d_o = (k * n) + d_o      # and denominator
            end
            
            rc[i * "|" * t] = c
            cf[i * "|" * t] = (rc[i * "|" * t] ) / d 

            rc_o[i * "|" * t] = c_o
            cf_o[i * "|" * t] = (rc_o[i * "|" * t] ) / d_o 
            
            ccn = i * "=" * f[i] * "|" * t * "=" * ts     # reconstruct cond class prob name
            ccn_o = i * "=" * f[i] * "|" * t * "=!" * ts  # and opposite
            
            println(rpad("    $ccn", 45), rpad("$c / $d", 15),  round(cf[i * "|" * t], digits=digits) )
            println(rpad("    $ccn_o", 45), rpad("$c_o / $d_o", 15),  round(cf_o[i * "|" * t], digits=digits), "\n" )

        end
        
        Pyx = tf * prod(values(cf))      # P(Y|X)
        Pnyx = tf_o * prod(values(cf_o)) # P(!Y|X)
        nd = Pyx + Pnyx                  # norm. denom.

        Pyx_norm = Pyx / nd
        Pnyx_norm = Pnyx / nd
        ARR = Pyx_norm/tf
        
        println("Probability of hypothesis\n")
        println("    P(Y|X): ", rpad(round(Pyx, digits=digits), 33), rpad("Normalized:", 15), round(Pyx_norm, digits=digits), "\n")
        
        println("Probability of alternative hypothesis\n")
        println("    P(!Y|X): ", rpad(round(Pnyx, digits=digits), 32), rpad("Normalized: ", 15),  round(Pnyx_norm, digits=digits), "\n")

        println("Absolute risk ratio [P(query) / baseline]\n")
        println(rpad("    P(Y|X) / P(Y):", 60), round(ARR, digits=digits))
        println(line)
        println()
        println()
        println(line)

        return Pyx_norm, Pnyx_norm, ARR, Pyx, Pnyx

    end

    
    Q1_Pyx_norm, Q1_Pnyx_norm, Q1_ARR, Q1_Pyx, Q1_Pnyx= getcounts(df, q1t, q1ts, q1f, k) 

    function boot_getcounts(dfin::DataFrame, target::String, targetstate::String, features::Dict, k::Int64, bootstraps=bootstraps)

        t = target; ts = targetstate; f=features; k=k;

        boot_Pyx_norm = []; boot_Pnyx_norm =  []; boot_ARR = []
        #boot_Pyx = []; boot_Pnyx = []
        
        for i in 1:bootstraps

            dfx = copy(dfin);    
            dfb = dfx[rand(1:nrow(dfx), size(dfx,1)), :]           

            rc   = Dict()    # raw counts
            cf   = Dict()    # corrected freq
            rc_o = Dict()    # raw counts !Y
            cf_o = Dict()    # corrected freq !Y

            rows = size(dfb,1)
            tc = sum(dfb[!,t] .== ts) # target counts given target state
            tf = tc / rows           # target state freq
            tc_o = sum(dfb[!,t] .!= ts)   # target counts opp
            tf_o = tc_o / rows           # target state freq opp

            for i in keys(f)         # class cond. probs
 
                fs = f[i]   
                c = sum(dfb[!, i] .== fs .&& dfb[!, t] .== ts)   # counts of Xi|Y
                d = sum(dfb[!, t] .== ts)                       # counts of Y where t=ts
                
                c_o = sum(dfb[!, i] .== fs .&& dfb[!, t] .!= ts)   # counts of Xi|Y
                d_o = sum(dfb[!, t] .!= ts)                       # counts of Y where t=ts
            
                if c == 0                       
                    n = length(unique(dfb[!, i]))
                    c = k + c            # class cond. prob correction for numerator
                    d = (k * n) + d      # and denominator
                end
            
                if c_o == 0
                    n = length(unique(dfb[!, i]))
                    c_o = k + c_o            # class cond. prob correction for numerator
                    d_o = (k * n) + d_o      # and denominator
                end
            
                rc[i * "|" * t] = c
                cf[i * "|" * t] = (rc[i * "|" * t] ) / d 
                
                rc_o[i * "|" * t] = c_o
                cf_o[i * "|" * t] = (rc_o[i * "|" * t] ) / d_o 
            
                ccn = i * "=" * f[i] * "|" * t * "=" * ts     # reconstruct cond class prob name
                ccn_o = i * "=" * f[i] * "|" * t * "=!" * ts  # and opposite

            end
        
            Pyx = tf * prod(values(cf))      # P(Y|X)
            Pnyx = tf_o * prod(values(cf_o)) # P(!Y|X)
            nd = Pyx + Pnyx                  # norm. denom.


            Pyx_norm = Pyx / nd
            Pnyx_norm = Pnyx / nd
            ARR = Pyx / tf
            push!(boot_Pyx_norm, Pyx_norm )
            push!(boot_Pnyx_norm, Pnyx_norm )
            push!(boot_ARR, ARR )

        end
        
        return boot_Pyx_norm, boot_Pnyx_norm, boot_ARR

    end
        
    if length(q2) > 0
        q2t = 0; q2ts = 0; q2f=Dict()  
        q2t, q2ts, q2f = getinput(q2)
        Q2_Pyx_norm, Q2_Pnyx_norm, Q2_ARR, Q2_Pyx, Q2_Pnyx = getcounts(df, q2t, q2ts, q2f, k)
  
        RRR = Q1_Pyx_norm / Q2_Pyx_norm
        println("Relative risk ratio [P(query1) vs. P(query2), normalized]\n")
        println(rpad("    Query1 P(Y|X) / Query2 P(Y|X)", 60),  round(RRR, digits=digits))
    end

    if bootstraps > 0

        boot_Q1_Pyx_norm=[]; boot_Q1_Pnyx_norm=[]; boot_Q1_ARR=[]
        boot_Q1_Pyx_norm, boot_Q1_Pnyx_norm, boot_Q1_ARR = boot_getcounts(df, q1t, q1ts, q1f, k, bootstraps)
#println(round.(boot_Q1_Pyx_norm, digits=4))
        boot_Q2_Pyx_norm=[]; boot_Q2_Pnyx_norm=[]; boot_Q2_ARR=[]
        boot_Q2_Pyx_norm, boot_Q2_Pnyx_norm, boot_Q2_ARR = boot_getcounts(df, q2t, q2ts, q2f, k, bootstraps)
#println(round.(boot_Q2_Pyx_norm, digits=4))
        boot_RRR = boot_Q1_Pyx_norm ./ boot_Q2_Pyx_norm 
        boot_RRR_med = median(boot_RRR)
        RRR05 = quantile(boot_RRR, 0.025)
        RRR95 = quantile(boot_RRR, 0.975)

        println(rpad("    Bootstrap median value:", 60), round(boot_RRR_med, digits=digits))
        println(rpad("    Bootstrap (CI95):", 55), "(", round(RRR05, digits=digits), ", ", round(RRR95, digits=digits), ")")
        println(line)
    else
        println(line)
    end

end
