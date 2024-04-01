"""
    cpq(data, "P(Y=y|X1=x1,X2=x2)")

Calculate an outcome probablity for the event Y=y given one or more conditional variables X1, X2, ... Xn from some body of evidence X.

Inputs: a csv file or dataframe and the probability query. The input data must have sample ids in column 1, and all columns must be labeled. Provide optional metatdata for samples in a separate file with sample ids in column 1.

Probabilities are estimated directly from the co-occurence of events in the data. The expression, P(Y=Yes|X1=Yes,X2=No), is queried as fraction of Y=Yes events in the subset of the data where X1 is Yes and X2 is No. Currently, only binary conditional variables are allowed. Laplace correction is use for zero count events.

    Keyword options:
        digits           Round to n digits                      [4]    
        outfile          Output file (appendable)               ""  
        outfile2         Simple output of ARR or RRR            ""   
        bootstraps       Number of replicates                   [100]  
        k                Laplace correction factor              [1.0]  
        mincounts        min target and conditional count flag  [1, 20]
        showids          Show target ids for query              false
        idinfo           Show additional info for target ids    ""
                         Input is csv file of metadata for ids 
        rrrdenom         Custom rel. risk ratio coditional      []
                         states for the denominator
                         (e.g., ["Yes", "No", "Yes"])
                         default: all conditional states negated
        confmeth         Distribution model for CI95        "empirical|t-dist"
        showhist         Display the ARR distribution            false 
        binomial         Binomial test of target counts          false
        fisher           Fisher test for ARR and RRR             false

Example query:  cpq(df, "P(Play=Yes|Forecast=Sunny,Wind=Yes)")
                The conditional probability of Play given Sunny and Wind.

Example query2: cpq(df, "P(Play=Yes, Forecast=Sunny|Wind=Yes)")
                The conditional probability of Play & Sunny given Wind.
                The Play and Forcast are merged into a single variable.

Notes: The Fisher and Binomial tests may be used initially as general guide for the interpretion of results, howevet, the underlying assumptions of variable independence is violated. Using the t-distribution is only recommend for queries with target size (e.g. < 30).

"""
function cpq(data::Union{String,DataFrame}, query::String; digits::Int64=4, k::Number=1, outfile::String="", bootstraps::Int64=0, mincounts::Array=[1,20], showids=false, outfile2::String="", rrrdenom::Array=[], confmeth="empirical", binomial::Bool=false, fisher::Bool=false, idinfo::String="", showhist::Bool=false)

    digits > 8 ? error("Significant digits maximum is 8 digits.") : "";
    line = "-"^72
    println(line)
    
    z = queryparser(query)

    newcols::AbstractVector{Symbol} = []
    push!(newcols, Symbol(split(z[1][1], "=")[1]))

    for i in z[2]
        x =  Symbol(split(i, "=")[1])
        push!(newcols, x)
    end
    
    if length(rrrdenom) > 0
        if length(rrrdenom) != length(z[2])
            error("\nrrrdenom entries must match the conditional variables.\n")
        end
        println("Using relative risk denominator states: $rrrdenom.")
    end
    

    if typeof(data) == DataFrame     #read data

        ids = data[!,1]
        nn = join(z[3], "")
        length(z[1]) > 1 ? push!(newcols, Symbol(nn)) : nothing
        
        if length(z[1]) > 1 && !in(nn, names(data)) 

            nnc = Symbol(nn)
            nd = Int64.(ones(size(data,1)))
            vu = string.(unique(data[!, Symbol(z[3][1])] ))
            
            if length(z[1]) == 2

                for row in eachindex(nd)
                    if ( data[row, Symbol(z[3][1])] == z[4][1] ) && ( data[row, Symbol(z[3][2])] == z[4][2] ) 
                        nd[row] = 2
                    end
                end
                
            elseif length(z[1]) == 3

                for row in eachindex(nd)
                    if ( data[row, Symbol(z[3][1])] == z[4][1] ) && ( data[row, Symbol(z[3][2])] == z[4][2] ) &&  ( data[row, Symbol(z[3][3])] == z[4][3] )
                        nd[row] = 2
                    end
                end

            else
            end

            nd = string.(nd)
            tv = String(z[4][1])
            
            if tv == vu[1]                
                nd = replace.(nd, "2" => vu[1], "1" => vu[2])
            else
                nd = replace.(nd, "2" => vu[2], "1" => vu[1])
            end                             
                 
            insertcols!(data, nn => nd)

        elseif length(z[4]) > 3 
            error("\n\nINFO: Multistate target evaluation limited to three targets!\n\n")   
        end

        df = select(data, newcols)

    elseif isfile(data)
        
        length(z[1]) > 1 ? error("\n\nDataframe input required for multitarget query\n\n.") : nothing
        
        global df = CSV.read(data, DataFrame, normalizenames=true, types=String, select=newcols);

        ids = df[!,1]

    else
        error("\nCannot find file or dataframe called $data.\n")
    end

    global dfg = df  #MUST set a global dataframe here for resampling with eval!!!
    global cc = 0
    global bb = -1

    sinfo=[]

    
    function runbnecpq(df, newcols, z)

        target = ""; tstate = "";

        for i in eachindex(z[1])

            tt = split(z[1][i], "=")
            target = String.(tt[1])
            tstate = String.(tt[2])

            if length(z[1]) > 1
                target = nn
                tstate = z[4][1]
            end

            if !in(tstate, df[!, target])
                tv = unique(df[!, target])
                if cc == 0
                    printstyled("\nWARNING: Target state $tstate not found in $target:\nObserved variables states: $tv.\n", color=:cyan)
                end
            end
        
        end

        features=""; fstates=""; q=""; q_opp=""

        for i in eachindex(z[2])      # z[2] has all conditionals in pairs

            ff = split(z[2][i], "=")
            feature = String.(ff[1])
            fstate = String.(ff[2])

            if !in(fstate, df[!, feature])
                fv = unique(df[!, feature])
                ft = typeof(fstate)
                printstyled("\nWARN: The feature state $fstate ($ft) was not found in $feature\nfor the query or a bootstrap query...\n", color=:yellow)
            end

            q = q * "df." * feature * ".==" * "\"" * fstate * "\"" * ","

            if length(rrrdenom) > 0          # custom array for rel risk denom.
                fstate_opp = [rrrdenom[i]]
            else
                fstate_opp = unique(df[!,feature][findall(df[!, feature] .!= fstate)])
            end

            if length(fstate_opp) == 1
                fstate_opp = fstate_opp[1]
            else
                error("\nMore than two states in selected feature\nor input inconsisent with features and states.\n")
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
        end
        
        targ_f = targ_c / all_c                    # target abs risk (baseline)
        t_feat_f = t_feat_c / cond_c               # targ cond freq 
        t_feat_f_fopp = t_feat_c_fopp / cond_c_opp # freq of target w/o conds
    
        # Use pseudocounts and simple Laplace correction (k=1 default)
        # for zero class problems. That is, given some feature/label,
        # if there are no target instances, 
        # correct zero probabablity to >0 and normalize by the number of
        # states of the target state.

        states = length(unique(df[!, target]))


        function laplace_correction(k, states, all_c, cond_c, cond_c_opp, targ_c, t_feat_c, t_feat_c_fopp)

            knorm =  k * states                     # laplace correction, norm by var states

            if t_feat_c == 0 || t_feat_c_fopp == 0  # lp correction of num and denom if
                t_c_adj = t_feat_c + k              # zero counts detected
                a_c_adj = cond_c + knorm        
                t_c_opp_adj = t_feat_c_fopp + k
                a_c_opp_adj = cond_c_opp + knorm
            else
                t_c_adj = t_feat_c                  # Don't correct non zero states
                a_c_adj = cond_c        
                t_c_opp_adj = t_feat_c_fopp 
                a_c_opp_adj = cond_c_opp
            end                

            lp_targ_f = targ_c / all_c                   # final frequencies
            lp_t_feat_f = t_c_adj / a_c_adj              # P(Y|X))
            lp_t_feat_f_fopp = t_c_opp_adj / a_c_opp_adj # P(Y|!X)

            lp_abs_risk_ratio = lp_t_feat_f / lp_targ_f
            lp_rel_risk_ratio =  lp_t_feat_f / lp_t_feat_f_fopp

            return (lp_targ_f, lp_t_feat_f, lp_t_feat_f_fopp, lp_abs_risk_ratio, lp_rel_risk_ratio, t_c_adj, a_c_adj, t_c_opp_adj, a_c_opp_adj)

        end

        lp_targ_f, lp_t_feat_f, lp_t_feat_f_fopp, lp_abs_risk_ratio, lp_rel_risk_ratio, t_c_adj, a_c_adj, t_c_opp_adj, a_c_opp_adj = laplace_correction(k, states, all_c, cond_c, cond_c_opp, targ_c, t_feat_c, t_feat_c_fopp)

        abs_risk_ratio = t_feat_f / targ_f
        rel_risk_ratio =  t_feat_f / t_feat_f_fopp

        fmt = join(["%.", digits, "f"], "")  # printf format string
        btpval = missing
        ftpvals = [missing, missing, missing]


        
        if cc == 0   #print the exact, non-bootstrapped counts and frequencies

            targ_f, t_feat_f, lp_t_feat_f,  t_feat_f_fopp, lp_t_feat_f_fopp = rpad.(printformat.(round.((targ_f, t_feat_f, lp_t_feat_f, t_feat_f_fopp, lp_t_feat_f_fopp ), digits=digits), fmt=fmt), 10)
            
            abs_risk_ratio, lp_abs_risk_ratio = rpad.(printformat.(round.((abs_risk_ratio,lp_abs_risk_ratio), digits=digits), fmt=fmt), 10)
            
            rel_risk_ratio, lp_rel_risk_ratio = rpad.(printformat.(round.((rel_risk_ratio,lp_rel_risk_ratio), digits=digits), fmt=fmt), 10)
            println(line)
            println("Counts:\n\tTarget (all)\t\t$targ_c\n\tTarget|Conditionals\t$t_feat_c\tConditionals only\t$all_c_cond\n\tTarget|!Conditionals\t$t_feat_c_fopp\t!Conditionals only\t$all_c_opp\n" )   
            println("Probabilities:\n  Absolute risk (baseline)\t$targ_f\n  Targets|Conditionals          $t_feat_f\t$lp_t_feat_f (adj)\n  Targets|!Conditionals         $t_feat_f_fopp\t$lp_t_feat_f_fopp (adj)\n")    
            println("Risk ratios:\n  ARR:\t$abs_risk_ratio\tARR (adj):\t$lp_abs_risk_ratio")
            println( "  RRR:\t$rel_risk_ratio \tRRR (adj):\t$lp_rel_risk_ratio")
            
            if (t_feat_c < mincounts[1] || all_c_cond < mincounts[2]) || (t_feat_c_fopp < mincounts[1] || all_c_opp < mincounts[2])
                println(line)
                error("\n\nInsufficient raw data counts using mincounts = $mincounts\nChange mincounts to run this query.\n\n")
                
            end

            if t_feat_c == 0
                printstyled("WARN: zero final counts observed. These results are only a prediction!\n", color=:yellow)
            end
            
            cc = cc + 1
            println(line)

            btpval = missing
            ftpvals = [missing, missing, missing, missing, missing, missing]

            if all_c_cond == 0  #Exact tests: Binomial and Fisher

                println("\n\nWARN: Zero conditional counts!\nSkipping exact tests for $query.\nAre the requested conditionals mutually exclusive?\n\n") 
                
            else

                dd = 1/10^digits

                if binomial == true
                    printstyled("-----------------------Binomial test of significance:-------------------\n", color=:cyan)
                    println(line)
                    bt = BinomialTest(t_feat_c, all_c_cond, parse.(Float64, targ_f))
                    btpval = round(pvalue(BinomialTest(t_feat_c, all_c_cond, parse.(Float64, targ_f))), digits=digits)
                    #bt95int = round.(confint(BinomialTest(t_feat_c, all_c_cond, parse.(Float64, targ_f))),  digits=digits)
                    btpval <= dd ? btpval = dd : btpval = btpval
                    
                    println(bt)
                    println(line)
                end
                
                if fisher == true
                
                    printstyled("---------------- Absolute risk ratio: Fisher evaluation ----------------\n", color=:green)
                    println(line)
                    
                    numall = all_c_cond - t_feat_c 
                    arrdenom = all_c - targ_c
                    
                    arrft = FisherExactTest(t_feat_c, targ_c, numall, arrdenom);

                    arrpvall = round(pvalue(FisherExactTest(t_feat_c, targ_c, numall, arrdenom), tail=:left ), digits=digits)
                    arrpvalr = round(pvalue(FisherExactTest(t_feat_c, targ_c, numall, arrdenom), tail=:right ), digits=digits)
                    arrpvalb = round(pvalue(FisherExactTest(t_feat_c, targ_c, numall, arrdenom), tail=:both ), digits=digits)

                    arrpvall <=  dd ? arrpvall = dd : arrpvall = arrpvall
                    arrpvalr <= dd ? arrpvalr = dd : arrpvalr = arrpvalr
                    arrpvalb <= dd ? arrpvalb = dd : arrpvalb = arrpvalb

                    println(arrft)
                    println("P-values:\n    Both-tails:   $arrpvalb\n    Left-tailed:  $arrpvall\n    Right-tailed: $arrpvalr\n\n")
                    println(line)
                    
                    #Fisher relative risk ratio
                    printstyled("------------ Relative risk ratio: Fishers Exact evaluation -------------\n", color=:cyan)
                    println(line)
                    numall = all_c_cond - t_feat_c
                    denomall = all_c_opp - t_feat_c_fopp
                    #println("FisherExactTest($t_feat_c, $t_feat_c_fopp, $numall, $denomall);")

                    ft = FisherExactTest(t_feat_c, t_feat_c_fopp, numall, denomall);

                    rrrpvall = round(pvalue(FisherExactTest(t_feat_c, t_feat_c_fopp, numall, denomall), tail=:left ), digits=digits)
                    rrrpvalr = round(pvalue(FisherExactTest(t_feat_c, t_feat_c_fopp, numall, denomall), tail=:right ), digits=digits)
                    rrrpvalb = round(pvalue(FisherExactTest(t_feat_c, t_feat_c_fopp, numall, denomall), tail=:both ), digits=digits)
                    
                    rrrpvall <= dd ? rrrpvall = dd : rrrpvall = rrrpvall
                    rrrpvalr <= dd ? rrrpvalr = dd : rrrpvalr = rrrpvalr
                    rrrpvalb <= dd ? rrrpvalb = dd : rrrpvalb = rrrpvalb

                    println(ft)
                    println("P-values:\n    Both-tails:   $rrrpvalb\n    Left-tailed:  $rrrpvall\n    Right-tailed: $rrrpvalr")
                    println(line)

                    ftpvals = [ arrpvalb, arrpvall, arrpvalr, rrrpvalb, rrrpvall,  rrrpvalr] 

                end
            end
            
            if showids == true 

                println("Observations driving the query probability:")
                println(tids)
                
                if length(idinfo) > 0

                    dfi = CSV.read(idinfo, DataFrame, types=String);
                    tids_a = split(tids, ",")
                    sinfo = DataFrame([name => [] for name in names(dfi)])

                    badids = setdiff(tids_a, data[:,1])
                        
                    for i in tids_a
                        info = dfi[dfi[!,1] .== i , :]
                        if size(info,1) > 0
                            sinfo = vcat(sinfo, info)
                        end                        
                    end
                    
                    println(line)
                    println("Information for identified observations:")
                    println(sinfo)

                    if length(badids) > 1
                        printstyled("WARN: Sample ids $badids were not found in the info file $idinfo.\n", color=:yellow)
                    end

                end
                println(line)
            end
        end
        
        return (([targ_f abs_risk_ratio lp_abs_risk_ratio rel_risk_ratio lp_rel_risk_ratio], [targ_c, all_c], [t_feat_c cond_c], [t_feat_c_fopp, all_c_opp], [t_c_adj, a_c_adj, t_c_opp_adj, a_c_opp_adj], [tids] ) ), btpval, ftpvals
        
    end

    rvals, btpval, ftpvals = runbnecpq(df, newcols, z)  # Main run

    if rvals[3][1] > 30 && confmeth == "t-dist"
        error("\n\nFinal target count > 30  for this probability query.\nPlease use confmeth=\"empirical\"\n\n")
    end
    
    
    function bnecpqbootstrap(dfg, bootstraps, rvals, btpval, ftpvals) # Input must be original df

        if confmeth == "t-dist" || confmeth == "empirical"
        else
            error("\n\nconfmeth must be \"t-dist\" or \"empirical\"\n\n")
        end
        
        HOUT = open("header.txt", "w")
        if length(rrrdenom) > 0
            println(HOUT, "##IMPORTANT: Denominator for relative risk calculations was set to $rrrdenom")
        end
        println(HOUT, "Query\tCount\tCondCount\tcpq_ARR_est\tcpq_ARR_CI95l\tcpq_ARR_CI95u\tcpq_RRR_est\tcpq_RRR_CI95l\tcpq_RRR_CI95u\tBinomial_pval\tARR_Fisher_2t_pval\tARR_Fisher_1t_left\tARR_Fisher_1t_right\tRRR_Fisher_2t_pval\tRRR_Fisher_1t_left\tRRR_Fisher_1t_right")
        close(HOUT)
        
        if bootstraps > 0
            
            if bootstraps < 100
                error("\n\nPlease use at least 100 bootstraps.\n\n")
            end
            
            B = Array{Float64, 2}(undef, bootstraps, 5)

            for i in 1:bootstraps
                df = dfg[rand(1:nrow(dfg), size(dfg,1)), :]
                b = runbnecpq(df, newcols, z)
                b = b[1]
                B[i,:] = b[1]              #2D array of resampled values
                cc = cc + 1
            end

            replace!(B, Inf=>NaN)  #set Inf to NaN and then filter invalid bootstraps
            
            B = sort(B, dims=1)
            B1 = filter(!isnan, B[:,1])
            B2 = filter(!isnan, B[:,2])
            B3 = filter(!isnan, B[:,3])
            B4 = filter(!isnan, B[:,4])
            B5 = filter(!isnan, B[:,5])

            function myconfint(data, bootstraps, digits)

                if length(data) < (bootstraps * 0.5)
                    e = -0.0; lc = -0.0; uc = -0.0
                   
                else

                    e = round(mean(data), digits=digits)

                    if confmeth == "t-dist"

                        #x = histogram(data)
                        #display(x)
                        #error("STOP")
                        N = rvals[3][1] # ν based on target count
                        N < 2 ? N = 2 : N = N # prevent ν domian error
                        
                        estd = std(data)
                        ci = confint(OneSampleTTest(e, estd, N), level=0.95)
                        lc = round.(ci[1], digits=digits)
                        uc = round.(ci[2], digits=digits)
                        lc < 0.0 ? lc = 0.0 : lc = lc
                    else
                        lc = quantile(data, 0.025)
                        uc = quantile(data, 0.975)
                    end

                    if showhist == true
                        plt = histogram(data, bins=50, xlabel="ARR score", ylabel="Observations ($bootstraps)", label="ARR estimates")
                        vline!([e], label="ARR average ($e)", width=1.0)
                        display(plt)
                    end
                    
                end
                
                return e, lc, uc                

            end

            TF, TFlc, TFuc = myconfint(B1, bootstraps, digits)   #get estimates and CI95s
            lp_ARR, lp_ARRlc, lp_ARRuc = myconfint(B3, bootstraps, digits)
            RRR, RRRlc, RRRuc = myconfint(B4, bootstraps, digits)
            lp_RRR, lp_RRRlc, lp_RRRuc = myconfint(B5, bootstraps, digits)

            ARR, ARRlc, ARRuc = myconfint(B2, bootstraps, digits)  #ARR last for showhist
            
            if cc == bootstraps + 1

                fmt = join(["%.", digits, "f"], "")
                
                TF, ARR, RRR, lp_ARR, lp_RRR = rpad.(printformat.(round.((TF, ARR, RRR, lp_ARR, lp_RRR ), digits=digits), fmt=fmt), 15)
                TFlc, TFuc, ARRlc, ARRuc, RRRlc, RRRuc = rpad.(printformat.(round.((TFlc, TFuc, ARRlc, ARRuc, RRRlc, RRRuc), digits=digits), fmt=fmt), 0)
                lp_ARRlc, lp_ARRuc, lp_RRRlc, lp_RRRuc = rpad.(printformat.(round.((lp_ARRlc, lp_ARRuc, lp_RRRlc, lp_RRRuc), digits=digits), fmt=fmt), 0)

                
                if bb > 0
                    println("WARN: $bb/$bootstraps bootstraps had zero counts but were corrected.")
                    println("-"^70)
                end
                
                println("Distribution average and CI95 for $bootstraps bootstraps")
                println("-"^70)
                println("Abs risk        $TF($TFlc, $TFuc)")
                println("ARR             $lp_ARR($lp_ARRlc, $lp_ARRuc)")
                println("RRR             $lp_RRR($lp_RRRlc, $lp_RRRuc)")
                println("-"^70)

                if (length(outfile) > 0) || (length(outfile2) > 0)

                    targraw     = rvals[2][1]
                    all_c       = rvals[2][2]
                    targcounts  = rvals[3][1]
                    condcounts  = rvals[3][2]
                    targopp     = rvals[4][1]
                    oppall      = rvals[4][2]
                    t_c_adj     = rvals[5][1]
                    a_c_adj     = rvals[5][2]
                    t_c_opp_adj = rvals[5][3]
                    a_c_opp_adj = rvals[5][4]

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

                    if condcounts < mincounts[2] || targcounts < mincounts[1]
                        badflag = "Insufficient_data"
                    end
                    
                    if targcounts == 0     # flag even if mincount set to 0
                        badflag = "Zero_raw_final_counts"
                    end
                    

                    for i in eachindex(clist)

                        if clist[i] >= mincounts[2]
                            vc[i] = 1
                        else
                            vc[i] = 0
                        end

                        if tlist[i] >= mincounts[1]
                            vt[i] = 1
                        else
                            vt[i] = 0
                        end

                    end

                    v1 = vc[1]; v2 = vc[2]; v3=vc[3]; v4 = vc[4]
                    mintarg1 = vt[1]; mintarg2 = vt[2]; mintarg3 = vt[3]; mintarg4 = vt[4]; 

                    targraw >= Int64(mincounts[1]) ? mintf = 1 : mintf = 0
                    targraw >= Int64(mincounts[2]) ? mintfc = 1  : mintfc = 0

                    (mintf == 0 || mintfc == 0) ? minsig = 0 : minsig = 1

                    
                    if length(outfile) > 0

                        if bootstraps < 1
                            error("Bootstrapping must be used for file output.")
                        end

                        # Binomial test and FisherTest vals
                        arrb = ftpvals[1]; arrl = ftpvals[2]; arrr = ftpvals[3]
                        rrrb = ftpvals[4]; rrrl = ftpvals[5]; rrrr = ftpvals[6]

                        OUT = open(outfile, "a")
                        println(OUT, "$query\t$t_c_adj\t$a_c_adj\t$lp_ARR\t$lp_ARRlc\t$lp_ARRuc\t$lp_RRR\t$lp_RRRlc\t$lp_RRRuc\t$btpval\t$arrb\t$arrl\t$arrr\t$rrrb\t$rrrl\t$rrrr")
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

    end            # end bnecpqbootstrap

    bnecpqbootstrap(df, bootstraps, rvals, btpval, ftpvals)
    
    fids = ""

    if showids
        fids = String.(split(rvals[6][1], r",")) 
    end


    arr = parse.(Float64, string.(strip(rvals[1][3]) ))
    rrr = parse.(Float64, string.(strip(rvals[1][5]) ))

    return arr, rrr, fids, sinfo, rvals[3][1]   #rval[3][1] is final counts
    
end
