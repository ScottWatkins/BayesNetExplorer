"""
    RRcalculator(dfp::DataFrame; target_state::String, target_state_freq::Float64, variable_table::DataFrame=vt, mincounts::Vector=[0,1], minp::Float=0.0, maxp::Float=1.0)

Calculate the relative and absolute risk ratios from a probability table generated with the bne iterate function. Process and sort all targets and target conditions. Users may apply filters to focus the results from a large probability space .

**Input variables, dfp, target_state_freq, and variable_table are from the output of the bne() run.**

All conditional variables in the input table must have exactly two states (typically 1/2 or 0/1 or true/false). The target variable may be multistate. All required inputs are generated from the bne() function. These include the main conditional probability table, the target state, the observed frequency of the target state, and the variable table. 

Options
-------
        minp            minimum probability to report         [0.0001]
        maxp            maximum probability to report         [0.9999]
        filterstates    remove queries if conditional         ""
                        variable states contain the string
        keep            keep queries if the conditional       ""
                        variable contain this string
        mincounts       min target and conditional counts     [1, 20]
        scols           array of col indices to sort output   [9,5]

    Example:

        dfr = RRcalculator(dfp, target_state="Mortality_true", target_state_freq=0.01, variable_table=vf)

Notes:
1) mincounts is set to [0,1] by default, but users may need to adjust these values (see discussion section in the manual). Combinations of mutually exclusive (e.g., hot one encoded) combinations that produce zero counts will be removed from the output.


"""
function RRcalculator(cpt::DataFrame; target_state::String, target_state_freq::Float64, variable_table::DataFrame, filterstates::String="", keep::String="", minp::Float64=0.0001, maxp::Float64=0.9999,  mincounts::Array{Int,1}=[0,1], scols::Array{Int,1}=[9,5] )

    if length(target_state_freq) < 1
        error("You must provide the baseline frequency of the target state observed in the whole dataset.\nThis is the population frequency over all samples.\n")
    end

    if length(unique(cpt.ConditionalVariables)) < 2
        error("\n\nRRcalculator iterates over queries with multiple conditional variables.\nPlease run bne with more than one conditional variable first.\nFor a single conditional variable test, please use bne, relrisk=true,\nto estimate the absolute and relative risk ratios.\n\n")
    end
    
    df_cpt = copy(cpt)

    if in(target_state, names(cpt)) == false
        opt1 = string(names(cpt)[2])
        opt2 = string(names(cpt)[3])
        error("\n\n$target_state not found.\nPossible values are $opt1 and $opt2.\n\n")
    end


    df_cpt = sort!( df_cpt[ (minp .≤ df_cpt[!,Symbol(target_state)] .≤ maxp) , :], Symbol(target_state), rev=true)
    opp_prob = Union{Float64,Missing}[]

    for i in 1:size(df_cpt, 1)

        cvars = df_cpt[i, 4]
        
        opp_state = replace(df_cpt[i, 5], "1" => "2", "2" => "1")

        df_cptx = df_cpt[ .&(df_cpt.ConditionalVariables .== cvars, df_cpt.ConditionalStates .== opp_state), :]

        if size(df_cptx,1) > 0
            y = df_cptx[ 1, Symbol(target_state) ]
            push!(opp_prob, y)
        else
            push!(opp_prob, missing)   
        end
        
    end

    insertcols!(df_cpt, 4, "P($target_state|!X)" => opp_prob)
    rr = df_cpt[!, Symbol(target_state)] ./ opp_prob
    insertcols!(df_cpt, 5, :RelRiskRatio => rr)

    abs_r = df_cpt[!, Symbol(target_state)] ./ target_state_freq
    insertcols!(df_cpt, 5, :AbsRiskRatio => abs_r)

    df_cpt = dropmissing(df_cpt)
    df_tt = variable_table

    ht = Dict()

    for i in 1:size(df_tt, 1) #lookup for conditional state names
        q = df_tt[i,1] * "_" *string(df_tt[i,2])
        ht[q] = string(df_tt[i,5])
    end

    vnames = []

    for i in 1:size(df_cpt,1)    #match conditional state with names

        cvd = []; csd = []; 
        cvd = split(df_cpt[i,7], ",")
        csd = split(df_cpt[i,8], ",")
        x = ""

        for j in eachindex(cvd)
            p = cvd[j] * "_" * csd[j]
            x = x * "," * string(get(ht, p, 0))
        end

        push!(vnames, x[2:end])

    end

    headdat = string.(split(replace(readline("BN.header"), "\"" => ""), " "))
    global dfc = CSV.read("BN.data", DataFrame, delim=" ", header=headdat)


    condcounts, tcounts = getcounts(df_cpt, dfc, target_state)

    insertcols!(df_cpt, 8, :ConditionalStateNames => vnames)
    insertcols!(df_cpt, 7, :CondCount => condcounts)
    insertcols!(df_cpt, 7, :Count => tcounts)
    np_probs = tcounts ./ condcounts
    insertcols!(df_cpt, 7, :NonPropProb => np_probs)
    
    sort!(df_cpt, (scols .+ 1) , rev=false)

    if length(filterstates) > 0
        df_cpt = df_cpt[ [!i for i in occursin.( Regex("$filterstates") , df_cpt.ConditionalStateNames) ] , :]
    end

    if length(keep) > 0
        df_cpt = df_cpt[ [i for i in occursin.( Regex("$keep") , df_cpt.ConditionalVariables) ] , :]
    end

    scs = size(df_cpt,1)
    
    df_cpt = df_cpt[df_cpt.CondCount .≥ mincounts[2], :]
    df_cpt = df_cpt[df_cpt.Count     .≥ mincounts[1], :]

    ecs = scs - size(df_cpt,1)
    println("Omitting $ecs rows based on mincounts = $mincounts ...")
    
    if names(df_cpt)[2] == target_state
        df_cpt = select(df_cpt, Not(3))
    else
        df_cpt = select(df_cpt, Not(2))
    end

    nn = join(["P(", target_state, "|X)"], "")
    rename!(df_cpt, Symbol(target_state) => nn)
    
    df_cpt[!,2:6] = round.(df_cpt[!, 2:6], digits=4)

    rename!(df_cpt, :ConditionalVariables  => :CondVariables)
    rename!(df_cpt, :ConditionalStateNames => :CondStateNames)
    rename!(df_cpt, :ConditionalStates     => :CondStates)

    return df_cpt
    
end
    
