"""
    bnemle(data, "P(A=a|B=b)")

Calculate the maximum likelihood probablity and CI95 for a target given some conditional variables.
Inputs are a full csv file or dataframe and the probability query.
Currently, only binary conditional variables are allowed. Laplace corrected risk ratios are provided (rounded at print time). Outfiles can be made for bootstrap runs.

    Keyword options:
        digits           Round to n digits                 [4]    
        outfile          Output file                        ""  
        outfile2         Single column output               ""   
        bootstraps       Number of replicates            [100]  
        k                Laplace correction factor         [1]  
                         (Set k=2 for Wilson midpoint)
        mincount         min target and cond count flag [1,20]
        showids          show target ids for query       false
"""
function bnemle(data::Union{String,DataFrame}, query::String; digits::Int64=4, k::Number=1, outfile::String="", bootstraps::Int64=0, mincount::Array=[1,20], showids=false, outfile2::String="")

    if occursin(r"^P\(((.+)=(.+))\|((.+)=(.+))\)$", query)
        println("-"^70)
        println("Query:$query")
    else
        error("\n\nBad conditional query: \"$query\"\nQuery format is \"P(A=X|B=Y,C=Z,...)\"\n")
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
        ids = data[!,1]
        df = select(data, newcols)
        global df = string.(df)
    elseif isfile(data)
        ids = data[!,1]
        global df = CSV.read(data, DataFrame, normalizenames=true, types=String, select=newcols)
    else
        error("\nCannot find file or dataframe called $data.\n")
    end

    global dfg = df  #MUST set a global dataframe here for resampling with eval!!!
    global cc = 0
    global bb = -1
    
    function runbnemle(df, newcols, z)

        target = ""; tstate = "";

        for i in eachindex(z[1])

            tt = split(z[1][i], "=")
            target = String.(tt[1])
            tstate = String.(tt[2])
        
            if !in(tstate, df[!, target])
                tv = unique(df[!, target])
                if cc == 0
                    printstyled("\nWARNING: Target state $tstate not found in $target:\nObserved variables states: $tv.\n", color=:cyan)
                end
            end
        
        end

        features=""; fstates=""; q=""; q_opp=""

        for i in eachindex(z[2])      #z[2] has all conditionals in pairs

            ff = split(z[2][i], "=")
            feature = String.(ff[1])
            fstate = String.(ff[2])

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

        global dfid = hcat(ids, df)

        oc = replace(q[8:end-5], "df" => "dfid")
        qid = "dfid[.&(" * oc  * "," * "dfid." * target * ".==\"" * tstate * "\"), 1]"
        dcall = Meta.parse(qid)
        tids = join(eval(dcall), ",")

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

        if t_feat_c == 0
            bb = bb + 1
            #println("WARN: Zero observations in conditional subset...")
        end
        
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

            knorm =  k * states         # laplace correction considering num target states
            t_adj = targ_c + k          # add k to numerator
            a_adj = all_c + knorm       # add knorm to denom.
            t_c_adj = t_feat_c + k  
            a_c_adj = cond_c + knorm
            t_c_opp_adj = t_feat_c_fopp + k
            a_c_opp_adj = cond_c_opp + knorm
            
            lp_targ_f = (targ_c + k) / (all_c + knorm)
            lp_t_feat_f = (t_feat_c + k) / (cond_c + knorm)
            lp_t_feat_f_fopp = (t_feat_c_fopp + k) / (cond_c_opp + knorm)
            lp_abs_risk_ratio = lp_t_feat_f / lp_targ_f
            lp_rel_risk_ratio =  lp_t_feat_f / lp_t_feat_f_fopp

            return (lp_targ_f, lp_t_feat_f, lp_t_feat_f_fopp, lp_abs_risk_ratio,          lp_rel_risk_ratio, t_adj, a_adj, t_c_adj, a_c_adj, t_c_opp_adj, a_c_opp_adj)

        end

        lp_targ_f, lp_t_feat_f, lp_t_feat_f_fopp, lp_abs_risk_ratio, lp_rel_risk_ratio, t_adj, a_adj, t_c_adj, a_c_adj, t_c_opp_adj, a_c_opp_adj = laplace_correction(k, states, all_c, cond_c, cond_c_opp, targ_c, t_feat_c, t_feat_c_fopp)

        abs_risk_ratio = t_feat_f / targ_f
        rel_risk_ratio =  t_feat_f / t_feat_f_fopp

        line = "-"^70

        fmt = join(["%.", digits, "f"], "") #printf format string for printformat function
        
        if cc == 0   #print only exact, non-bootstrap values

            targ_f, t_feat_f, lp_targ_f, lp_t_feat_f,  t_feat_f_fopp, lp_t_feat_f_fopp = rpad.(printformat.(round.((targ_f, t_feat_f, lp_targ_f, lp_t_feat_f, t_feat_f_fopp, lp_t_feat_f_fopp ), digits=digits), fmt=fmt), 10)
            
            abs_risk_ratio, lp_abs_risk_ratio = rpad.(printformat.(round.((abs_risk_ratio,lp_abs_risk_ratio), digits=digits), fmt=fmt), 10)
            
            rel_risk_ratio, lp_rel_risk_ratio = rpad.(printformat.(round.((rel_risk_ratio,lp_rel_risk_ratio), digits=digits), fmt=fmt), 10)
            println(line)
            println("Counts:\n\tTarget (all)\t\t$targ_c\n\tTarget|Conditionals\t$t_feat_c\tConditionals only\t$all_c_cond\n\tTarget|!Conditionals\t$t_feat_c_fopp\t!Conditionals only\t$all_c_opp\n" )   
            println("Frequencies:\n  Absolute risk (baseline)\t$targ_f\t$lp_targ_f (adj)\n  Targets|Conditionals          $t_feat_f\t$lp_t_feat_f (adj)\n  Targets|!Conditionals         $t_feat_f_fopp\t$lp_t_feat_f_fopp (adj)\n")    
            println("Risk ratios:\n  ARR:\t$abs_risk_ratio\tARR (adj):\t$lp_abs_risk_ratio")
            println( "  RRR:\t$rel_risk_ratio \tRRR (adj):\t$lp_rel_risk_ratio")

            if all_c_cond < mincount[2] || targ_c < mincount[1]
                println(line)
                println("WARN: Insufficient raw data for this query (using mincount = $mincount)")
            end



            cc = cc + 1
            println(line)
            
            if showids == true 
                println("Observations driving the probabilities listed above:")
                println(tids)
                println(line)
            end
        end

        return(([targ_f abs_risk_ratio lp_abs_risk_ratio rel_risk_ratio lp_rel_risk_ratio], [targ_c, all_c], [t_feat_c cond_c], [t_feat_c_fopp, all_c_opp], [t_adj, a_adj, t_c_adj, a_c_adj, t_c_opp_adj, a_c_opp_adj], [tids] )  )


    end

    rvals = runbnemle(df, newcols, z)

    
    function bnemlebootstrap(dfg, bootstraps, rvals) # Input is original df

        HOUT = open("header.txt", "w")
        println(HOUT, "Query\tEstimate\tDistMean\tCI95l\tCI95u\tCondCount\tFinalCount\tTargCountOK\tCondCountOK\tSignificant\tQueryOK")
        close(HOUT)

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

                if bb > 0
                    println("WARN: $bb/$bootstraps bootstraps had zero count states (laplace corrected).")
                    println("-"^70)
                end
                
                println("Distribution average and CI95 for $bootstraps bootstraps")
                println("-"^70)
                println("Abs risk        $TF($TFlc, $TFuc)")
                println("ARR             $ARR($ARRlc, $ARRuc)\nARR (adj)       $lp_ARR($lp_ARRlc, $lp_ARRuc)")
                println("RRR             $RRR($RRRlc, $RRRuc)\nRRR (adj)       $lp_RRR($lp_RRRlc, $lp_RRRuc)")
                println("-"^70)
                
                if (length(outfile) > 0) || (length(outfile2) > 0)
                    
                    targraw     = rvals[2][1]
                    all_c       = rvals[2][2]
                    targcounts  = rvals[3][1]
                    condcounts  = rvals[3][2]
                    targopp     = rvals[4][1]
                    oppall      = rvals[4][2]
                    t_adj       = rvals[5][1]
                    a_adj       = rvals[5][2]
                    t_c_adj     = rvals[5][3]
                    a_c_adj     = rvals[5][4]
                    t_c_opp_adj = rvals[5][5]
                    a_c_opp_adj = rvals[5][6]

                    sig = Vector{Int64}(undef, 4) #generate a different than 1 flag
                    svals_lc = parse.(Float64, [ARRlc, RRRlc, lp_ARRlc, lp_RRRlc])
                    svals_uc = parse.(Float64, [ARRuc, RRRuc, lp_ARRuc, lp_RRRuc])

                    for i in eachindex(svals_lc)
                        if svals_lc[i] <= 1.0 < svals_uc[i]
                            sig[i] = 0
                        else
                            sig[i] = 1
                        end
                    end
                    
                    s1 = sig[1]; s2 = sig[2]; s3=sig[3]; s4 = sig[4]

                    vc = Vector{Int64}(undef, 4) #generate valid count flag
                    clist = [condcounts, a_c_adj, oppall, a_c_opp_adj]
                    vt = Vector{Int64}(undef, 4) #generate valid targ count flag
                    tlist = [targcounts, t_c_adj, targopp, t_c_opp_adj]

                    badflag = "OK"
                    if condcounts < mincount[2] || targcounts < mincount[1]
                        badflag = "Insufficient_data"
                    end
                    
                    for i in eachindex(clist)

                        if clist[i] >= mincount[2]
                            vc[i] = 1
                        else
                            vc[i] = 0
                        end

                        if tlist[i] >= mincount[1]
                            vt[i] = 1
                        else
                            vt[i] = 0
                        end

                    end


                    v1 = vc[1]; v2 = vc[2]; v3=vc[3]; v4 = vc[4]
                    mintarg1 = vt[1]; mintarg2 = vt[2]; mintarg3 = vt[3]; mintarg4 = vt[4]; 

                    targraw >= Int64(mincount[1]) ? mintf = 1 : mintf = 0
                    targraw >= Int64(mincount[2]) ? mintfc = 1  : mintfc = 0

                    (mintf == 0 || mintfc == 0) ? minsig = 0 : minsig = 1
                    if length(outfile) > 0
                        OUT = open(outfile, "a")
                        println(OUT, "$query\tAR\t$TF\t$TFlc\t$TFuc\t$all_c\t$targraw\t$mintf\t$mintfc\t$minsig\t$badflag")
                        println(OUT, "$query\tARR\t$ARR\t$ARRlc\t$ARRuc\t$condcounts\t$targcounts\t$mintarg1\t$v1\t$s1\t$badflag")
                        println(OUT, "$query\tARR (adj)\t$lp_ARR\t$lp_ARRlc\t$lp_ARRuc\t$a_c_adj\t$t_c_adj\t$mintarg2\t$v2\t$s2\t$badflag")
                        println(OUT, "$query\tRRR\t$RRR\t$RRRlc\t$RRRuc\t$oppall\t$targopp\t$mintarg3\t$v3\t$s3\t$badflag")
                        println(OUT, "$query\tRRR (adj)\t$lp_RRR\t$lp_RRRlc\t$lp_RRRuc\t$a_c_opp_adj\t$t_c_opp_adj\t$mintarg4\t$v4\t$s4\t$badflag")
                        close(OUT)
                    end

                    if length(outfile2) > 0
                        OUT2 = open(outfile2, "a")

                        if occursin(r"arr"i, outfile2)                                                      
                            println(OUT2, query, ",", ARR)
                        elseif occursin(r"rrr"i, outfile2)
                            println(OUT2, query, ",", RRR)
                        else
                            error("\n\nPlease include ARR or RRR in the outfile2 name.\n\n")
                        end

                        close(OUT2)
                    end
                    
                end
                
            end
        end
    end            # end bnemlebootstrap

    bnemlebootstrap(df, bootstraps, rvals)

    if showids
        fids = String.(split(rvals[6][1], r",")) 
        return fids
    end

end
