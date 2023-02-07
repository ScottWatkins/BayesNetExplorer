"""
    bnemle(data, "P(A=a|B=b)")

Calculate the maximum likelihood probablity and CI95 for a target given some conditional variables.
Inputs are a full csv file or dataframe and the probability query.
Currently, only binary conditional variables are allowed. Laplace corrected risk ratios are provided.

    Keyword options:
        digits           Round to n digits            [4]    
        outfile          Output file                   ""  
        bootstraps       Number of replicates       [100]  
        k                Laplace correction factor    [1]  
"""
function bnemle(data::Union{String,DataFrame}, query::String; digits::Int64=4, k::Number=1, outfile::String="", bootstraps::Int64=0)

    if occursin(r"^P\(((.+)=(.+))\|((.+)=(.+))\)$", query)
        println("-"^70)
        println("Query:\n$query")
    else
        error("\n\nBad query format: $query\n\n")
    end
    digits > 8 ? error("Significant digits should be less than 8.") :
        
    z = queryparser(query)

    newcols::AbstractVector{Symbol} = []
    push!(newcols, Symbol(split(z[1][1], "=")[1]))

    for i in z[2]
        x =  Symbol(split(i, "=")[1])
        push!(newcols, x)
    end

    if typeof(data) == DataFrame     #read data
        df = select(data, newcols)
        global df = string.(df)
    elseif isfile(data)
        global df = CSV.read(data, DataFrame, normalizenames=true, types=String, select=newcols)
    else
        error("\nCannot find file or dataframe called $data.\n")
    end

    global dfg = df  #MUST set a global dataframe here for resampling with eval!!!
    global cc = 0
    
    function runbnemle(df, newcols, z)

        target = ""; tstate = "";

        for i in eachindex(z[1])

            tt = split(z[1][i], "=")
            target = String.(tt[1])
            tstate = String.(tt[2])
        
            if !in(tstate, df[!, target])
                tv = unique(df[!, target])
                printstyled("\nWARNING: Target state $tstate not found in $target:\nVariables states: $v.\n", color=:cyan)
            end
        
        end

        features=""; fstates=""; q=""; q_opp=""

        for i in eachindex(z[2])      #z[2] has all conditionals in pairs

            ff = split(z[2][i], "=")
            feature = String.(ff[1])
            fstate = String.(ff[2])
#println("===>", feature, " is ", fstate)
            if !in(fstate, df[!, feature])
                fv = unique(df[!, feature])
                ft = typeof(fstate)
                printstyled("\nWARNING: The feature state $fstate ($ft) was not found in $feature:\nObserved variables states were: $fv.\nThe conditional query is not valid, but an adjusted value may
still be reported when this warning comes from a resampled dataset.\n", color=:yellow)
            end

            q = q * "df." * feature * ".==" * "\"" * fstate * "\"" * ","

            fstate_opp = unique(df[!,feature][findall(df[!, feature] .!= fstate)])

            if length(fstate_opp) == 1
                fstate_opp = fstate_opp[1]
            else
                error("More than two states in selected feature.")
            end

            q_opp = q_opp * "df." * feature * ".==" * "\"" * fstate_opp * "\"" * ","

        end

        dfc = DataFrame(); dfc_opp = DataFrame();
        
        q = "df[ .&(" * q[1:end-1] * "), :]"

        dcall = Meta.parse(q)
        dfc = eval(dcall)

        q_opp = "df[ .&(" * q_opp[1:end-1] * "), :]"

        dcallopp = Meta.parse(q_opp)
        dfc_opp = eval(dcallopp)
    
        all_c = size(df,1)
        vars = size(df,2)
        all_c_cond = size(dfc,1) 
        all_c_opp = size(dfc_opp,1)
        all_f_opp = all_c_opp/all_c

        targ_c = sum(df[!, target] .== tstate)    # total targ count all
        t_feat_c = sum(dfc[!, target] .== tstate) # targ in cond subset
        t_feat_c_fopp = sum(dfc_opp[!, target] .== tstate) # targ count w/o conds
        cond_c = size(dfc,1)                      # cond count no targ
        cond_c_opp = size(dfc_opp,1)              # cond count opp no targ

        targ_f = targ_c / all_c                   # target abs risk (baseline)
        t_feat_f = t_feat_c / cond_c              # targ cond freq 
        t_feat_f_fopp = t_feat_c_fopp / cond_c_opp  # freq of target w/o conds
    
        # With pseudocounts and simple Laplace correction (k=1 default)
        # for the zero class problem. That is, given some feature/label
        # the are no target instances for at least one state
        # so zero probabablity corrected to >0 and normalized by the number of
        # states of the target state.
        # https://classes.soe.ucsc.edu/cmps140/Winter17/slides/3.pdf

        states = length(unique(df[!, target]))

        function laplace_correction(k, states, all_c, cond_c, cond_c_opp, targ_c, t_feat_c, t_feat_c_fopp)
            knorm = states * k
            lp_targ_f = (targ_c + k) / (all_c + knorm)
            lp_t_feat_f = (t_feat_c + k) / (cond_c + knorm)
            lp_t_feat_f_fopp = (t_feat_c_fopp + k) / (cond_c_opp + knorm)
            lp_abs_risk_ratio = lp_t_feat_f / lp_targ_f
            lp_rel_risk_ratio =  lp_t_feat_f / lp_t_feat_f_fopp
            return (lp_targ_f, lp_t_feat_f, lp_t_feat_f_fopp, lp_abs_risk_ratio, lp_rel_risk_ratio)
        end

        lp_targ_f, lp_t_feat_f, lp_t_feat_f_fopp, lp_abs_risk_ratio, lp_rel_risk_ratio = laplace_correction(k, states, all_c, cond_c, cond_c_opp, targ_c, t_feat_c, t_feat_c_fopp)

        abs_risk_ratio = t_feat_f / targ_f
        rel_risk_ratio =  t_feat_f / t_feat_f_fopp

        line = "-"^70

        fmt = join(["%.", digits, "f"], "") #printf format string for printformat function
        
        if cc == 0   #print only exact, non-bootstrap values

            targ_f, t_feat_f, t_feat_f_fopp = rpad.(printformat.(round.((targ_f, t_feat_f, t_feat_f_fopp ), digits=digits), fmt=fmt), 10)
            abs_risk_ratio, lp_abs_risk_ratio = rpad.(printformat.(round.((abs_risk_ratio,lp_abs_risk_ratio), digits=digits), fmt=fmt), 10)
            rel_risk_ratio, lp_rel_risk_ratio = rpad.(printformat.(round.((rel_risk_ratio,lp_rel_risk_ratio), digits=digits), fmt=fmt), 10)
            println(line)
            println("Counts:\n\tTarget (all)\t\t$targ_c\n\tTarget|Conditionals\t$t_feat_c\tConditionals only\t$all_c_cond\n\tTarget|!Conditionals\t$t_feat_c_fopp\t!Conditionals only\t$all_c_opp\n" )   
            println("Frequencies:\n  Absolute risk (baseline)\t$targ_f\n  Targets|Conditionals          $t_feat_f\n  Targets|!Conditionals         $t_feat_f_fopp\n")    
            println("Risk ratios:\n  ARR:\t$abs_risk_ratio\tARR (adj):\t$lp_abs_risk_ratio")
            println( "  RRR:\t$rel_risk_ratio \tRRR (adj):\t$lp_rel_risk_ratio")
            cc = cc + 1
            println(line)

            #if length(outfile) > 0
            #    OUT = open(outfile, "a")
            #    o = join([query, targ_f, abs_risk_ratio, lp_abs_risk_ratio, rel_risk_ratio, lp_rel_risk_ratio], ",")
            #    println(OUT, o)
            #    close(OUT)
            #end

        end
        
        return(([targ_f abs_risk_ratio lp_abs_risk_ratio rel_risk_ratio lp_rel_risk_ratio], [t_feat_c cond_c]))

    end

    rvals = runbnemle(df, newcols, z)

    function bnemlebootstrap(dfg, bootstraps, rvals) # Input is original df
        
        if bootstraps > 0

            B = Array{Float64, 2}(undef, bootstraps, 5)
            
            for i in 1:bootstraps
                df = dfg[rand(1:nrow(dfg), size(dfg,1)), :]
                b = runbnemle(df, newcols, z)
                B[i,:] = b[1]              #2D array of rsampled values
                cc = cc + 1
            end
            
            lc = Int64(ceil(bootstraps * 0.025))    #CI95 cuts
            uc = Int64(floor(bootstraps * 0.975))
            B = sort(B, dims=1)
            
            TF = mean(B[:,1]);         TFlc = B[lc,1];     TFuc = B[uc,1];
            ARR = mean(B[:,2]);       ARRlc = B[lc,2];    ARRuc = B[uc,2];
            lp_ARR = mean(B[:,3]); lp_ARRlc = B[lc,3]; lp_ARRuc = B[uc,3];
            RRR = mean(B[:,4]);       RRRlc = B[lc,4];    RRRuc = B[uc,4];
            lp_RRR = mean(B[:,5]); lp_RRRlc = B[lc,5]; lp_RRRuc = B[uc,5];

            if cc == bootstraps + 1

                fmt = join(["%.", digits, "f"], "")
                
                TF, ARR, RRR, lp_ARR, lp_RRR = rpad.(printformat.(round.((TF, ARR, RRR, lp_ARR, lp_RRR ), digits=digits), fmt=fmt), 15)
                TFlc, TFuc, ARRlc, ARRuc, RRRlc, RRRuc = rpad.(printformat.(round.((TFlc, TFuc, ARRlc, ARRuc, RRRlc, RRRuc), digits=digits), fmt=fmt), 0)
                lp_ARRlc, lp_ARRuc, lp_RRRlc, lp_RRRuc = rpad.(printformat.(round.((lp_ARRlc, lp_ARRuc, lp_RRRlc, lp_RRRuc), digits=digits), fmt=fmt), 0)

                println("Distribution average and CI95 for $bootstraps bootstraps")
                println("-"^70)
                println("Abs risk        $TF($TFlc, $TFuc)")
                println("ARR             $ARR($ARRlc, $ARRuc)\nARR (adj)       $lp_ARR($lp_ARRlc, $lp_ARRuc)")
                println("RRR             $RRR($RRRlc, $RRRuc)\nRRR (adj)       $lp_RRR($lp_RRRlc, $lp_RRRuc)")
                println("-"^70)
                
                if length(outfile) > 0
                    targcounts = rvals[2][1]
                    condcounts = rvals[2][2]
                    
                    OUT = open(outfile, "a")
                    #println(OUT, "$query\tAR\t\t$TF\t($TFlc, $TFuc)")
                    #println(OUT, "$query\tARR\t\t$ARR\t($ARRlc, $ARRuc")
                    println(OUT, "$query\tARR (adj)\t$lp_ARR\t($lp_ARRlc, $lp_ARRuc)\t$condcounts\t$targcounts")
                    #println(OUT, "$query\tRRR\t\t$RRR\t($RRRlc, $RRRuc)
                    #println("$query\tRRR (adj)\t$lp_RRR\t($lp_RRRlc, $lp_RRRuc)")
                    close(OUT)
                end
                
            end
        end
    end

    bnemlebootstrap(df, bootstraps, rvals)
    
end
