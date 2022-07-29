"""
    RRcalculator(cpt::DataFrame; target_state::String, target_state_freq::Float64, translation_table::DataFrame=tt, mincondcount::Int=20, minp::Float=0.01, maxp=0.99)

Calculate the relative and absolute risk ratios from a probability table generated with the bne function. Process and sort all targets and target conditions. Users can apply filters to help remove low confidence conditional probabilities. All variables in the input table must have exactly two states (e.g. 0/1 or true/false). 

All required inputs are generated from the bne() function. These include the main conditional probability table, the target state, the observed frequency of the target state, and the variable translation table. 

Note: mincondcount is set to 20 by default and may need to be set to a lower value if a small data set is analyzed (see discussion section in the manual).

Options:

minp             minimum probability to report    [0.01]
maxp             maximum probability to report    [0.99]
filterstates           remove any row with any conditional variables sets containing this value    []
                 Typically Used to remove falses and see only states where the conditional is true
mincondcount         report results with at least this many samples in the *conditional* subset    [20]

Example:

RRcalculator(cpt, target\\_state="Mortality\\_true", target\\_state\\_freq=0.001, translation\\_table=dft, filterstates="", keep="",  mincondcount=20, minp=0.01, maxp=0.99 )



"""
function RRcalculator(cpt::DataFrame; target_state::String, target_state_freq::Float64, translation_table::DataFrame, filterstates::String="", keep::String="", minp::Float64=0.01, maxp::Float64=0.99,  mincondcount::Int=1)

    if length(target_state_freq) < 1
        error("You must provide the baseline frequency of the target state observed in the whole dataset.\nThis is the population frequency over all samples.\n")
    end

    if size(cpt, 1) < 3
        error("\n\nRRcalculator iterates over two or more conditional traits.\nPlease add additional traits, or use bne to estimate\nthe relative risk for a single conditional feature.\n\n")
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

    insertcols!(df_cpt, 4, "Prob(X)|(~Y)" => opp_prob)
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

    nn = join(["P(", target_state, ")|Y"], "")
    rename!(df_cpt, Symbol(target_state) => nn)
    
    df_cpt[!,2:6] = round.(df_cpt[!, 2:6], digits=4)

    return df_cpt
    
end
    
