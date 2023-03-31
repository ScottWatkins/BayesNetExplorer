"""
    RRcalculator(cpt::DataFrame; target_state::String, target_state_freq::Float64, translation_table::DataFrame=tt, mincondcount::Int=20, minp::Float=0.0, maxp::Float=1.0)

Calculate the relative and absolute risk ratios from a conditional probability table generated with the bne function. Process and sort all targets and target conditions. Users may apply filters to help remove low confidence conditional probabilities. All conditional variables in the input table must have exactly two states (e.g. 0/1 or true/false). The target variable can be multistate.

All required inputs are generated from the bne() function. These include the main conditional probability table, the target state, the observed frequency of the target state, and the variable translation table. 

Notes:
1) mincondcount is set to 20 by default and may need to be set to a lower value if the data set is small or the target variable state has a low frequency (see discussion section in the manual).
2) Combinations of  mutually exclusive (e.g., hot one encoded) combinations produce that produce zero counts will be omitted from the output.

    Options:

        minp            minimum probability to report    [0.00001]
        maxp            maximum probability to report    [0.99999]
        filterstates    remove conditional variables     ""
        keep            keep conditional variable        ""
        mincondcount    min conditional subset count     [20]

    Example:

        RRcalculator(cpt, target_state="Mortality_true", target_state_freq=0.01, translation_table=tt)

"""
function RRcalculator(cpt::DataFrame; target_state::String, target_state_freq::Float64, translation_table::DataFrame, filterstates::String="", keep::String="", minp::Float64=0.00001, maxp::Float64=0.99999,  mincondcount::Int=1)

    if length(target_state_freq) < 1
        error("You must provide the baseline frequency of the target state observed in the whole dataset.\nThis is the population frequency over all samples.\n")
    end

    if size(cpt, 1) < 2
        error("\n\nRRcalculator iterates over two or more target queries.\nPlease add additional queries, or use bne(... relrisk=true) to estimate\nthe relative risk for a single target query.\n\n")
    end
    
    df_cpt = copy(cpt)

    if in(target_state, names(cpt)) == false
        error("\n\n$target_state not found. Target_state format is feature_state (e.g. Arrhythmia_true or Cancer_false)\n\n")
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

    insertcols!(df_cpt, 4, "P($target_state|!Y)" => opp_prob)
    rr = df_cpt[!, Symbol(target_state)] ./ opp_prob
    insertcols!(df_cpt, 5, :RelRiskRatio => rr)

    abs_r = df_cpt[!, Symbol(target_state)] ./ target_state_freq
    insertcols!(df_cpt, 5, :AbsRiskRatio => abs_r)

    df_cpt = dropmissing(df_cpt)
    df_tt = translation_table

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
    
    sort!(df_cpt, [:RelRiskRatio, :CondCount], rev=true)

    if length(filterstates) > 0
        df_cpt = df_cpt[ [!i for i in occursin.( Regex("$filterstates") , df_cpt.ConditionalStateNames) ] , :]
    end

    if length(keep) > 0
        df_cpt = df_cpt[ [i for i in occursin.( Regex("$keep") , df_cpt.ConditionalVariables) ] , :]
    end

    df_cpt = df_cpt[df_cpt.CondCount .≥ mincondcount, :]

    if names(df_cpt)[2] == target_state
        df_cpt = select(df_cpt, Not(3))
    else
        df_cpt = select(df_cpt, Not(2))
    end

    nn = join(["P(", target_state, "|Y)"], "")
    rename!(df_cpt, Symbol(target_state) => nn)
    
    df_cpt[!,2:6] = round.(df_cpt[!, 2:6], digits=4)

    rename!(df_cpt, :ConditionalVariables  => :CondVariables)
    rename!(df_cpt, :ConditionalStateNames => :CondStateNames)
    rename!(df_cpt, :ConditionalStates     => :CondStates)

    return df_cpt
    
end
    
