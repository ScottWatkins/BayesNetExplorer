"""
    getRRtablerows(RRtable::DataFrame; targets::Array=[], condvars::Array=[], condstates::Array=[])

Retrieve conditinal query results from a relative risk table (RRtable) created by iterate="all" or "nodes".

this function simplies examining relative risk states when many variables are iterated using bne iteration.

Example: getRRtablerows(dfa, targets=["MORT","EFF"], condvars=["MUTATION","VENT"], condstates=["Yes","Yes])

"""
function getRRtablerows(RRtable::DataFrame; targets::Array=[], condvars::Array=[], condstates::Array=[])

    dfa = RRtable
    dfz = DataFrame()
    a = condvars[1];
    cs1 = condstates[1]

    for s in targets

        if length(condvars) == 1
            dfz = vcat(dfz, dfa[.&(dfa.Target .== s, .&(occursin.(Regex(a), dfa.CondVariables)), dfa.CondStates .== cs1), :])
        end

        if length(condvars) == 2
            b = condvars[2]; cs2 = condstates[2]; cs12 = join([cs1, cs2], ",");
            dfz = vcat(dfz, dfa[.&(dfa.Target .== s, .&(occursin.(Regex(a), dfa.CondVariables)), dfa.CondStates .== cs1), :])
            dfz = vcat(dfz, dfa[.&(dfa.Target .== s, .&(occursin.(Regex(b), dfa.CondVariables)), dfa.CondStates .== cs2), :])
            dfz = vcat(dfz, dfa[.&(dfa.Target .== s, .&(occursin.(Regex(a), dfa.CondVariables), occursin.(Regex(b), dfa.CondVariables)), dfa.CondStates .== cs12), :])
        end

        if length(condvars) == 3
            b = condvars[2]; c = condvars[3];
            cs2 = condstates[2]; cs3 = condstates[3];
            cs12 = join([cs1, cs2], ","); cs13 = join([cs1, cs3], ","); cs23 = join([cs2, cs3], ",");  cs123 = join([cs1, cs2, cs3], ",");
            dfz = vcat(dfz, dfa[.&(dfa.Target .== s, .&(occursin.(Regex(a), dfa.CondVariables)), dfa.CondStates .== cs1), :])
            dfz = vcat(dfz, dfa[.&(dfa.Target .== s, .&(occursin.(Regex(b), dfa.CondVariables)), dfa.CondStates .== cs2), :])
            dfz = vcat(dfz, dfa[.&(dfa.Target .== s, .&(occursin.(Regex(c), dfa.CondVariables)), dfa.CondStates .== cs3), :])
            dfz = vcat(dfz, dfa[.&(dfa.Target .== s, .&(occursin.(Regex(a), dfa.CondVariables), occursin.(Regex(b), dfa.CondVariables)), dfa.CondStates .== cs12), :])
            dfz = vcat(dfz, dfa[.&(dfa.Target .== s, .&(occursin.(Regex(a), dfa.CondVariables), occursin.(Regex(c), dfa.CondVariables)), dfa.CondStates .== cs13), :])
            dfz = vcat(dfz, dfa[.&(dfa.Target .== s, .&(occursin.(Regex(b), dfa.CondVariables), occursin.(Regex(c), dfa.CondVariables)), dfa.CondStates .== cs23), :])
            dfz = vcat(dfz, dfa[.&(dfa.Target .== s, .&(occursin.(Regex(a), dfa.CondVariables), occursin.(Regex(b), dfa.CondVariables), occursin.(Regex(c), dfa.CondVariables)  ), dfa.CondStates .== cs123), :])
        end

        if length(condvars) == 4
            b = condvars[2]; c = condvars[3]; d = condvars[4];
            cs2 = condstates[2]; cs3 = condstates[3]; cs4 = condstates[4];
            cs12 = join([cs1, cs2], ","); cs13 = join([cs1, cs3], ","); cs14 = join([cs1, cs4], ",");
            cs23 = join([cs2, cs3], ","); cs24 = join([cs2, cs4], ","); cs34 = join([cs3, cs4], ",");
            cs123 = join([cs1, cs2, cs3], ","); cs124 = join([cs1, cs2, cs4], ","); cs134 = join([cs1, cs3, cs4], ",");
            cs234 = join([cs2, cs3, cs4], ",");; cs1234 = join([cs1, cs2, cs3, cs4], ",");

            dfz = vcat(dfz, dfa[.&(dfa.Target .== s, .&(occursin.(Regex(a), dfa.CondVariables)), dfa.CondStates .== cs1), :])
            dfz = vcat(dfz, dfa[.&(dfa.Target .== s, .&(occursin.(Regex(b), dfa.CondVariables)), dfa.CondStates .== cs2), :])
            dfz = vcat(dfz, dfa[.&(dfa.Target .== s, .&(occursin.(Regex(c), dfa.CondVariables)), dfa.CondStates .== cs4), :])
            dfz = vcat(dfz, dfa[.&(dfa.Target .== s, .&(occursin.(Regex(d), dfa.CondVariables)), dfa.CondStates .== cs4), :])

            dfz = vcat(dfz, dfa[.&(dfa.Target .== s, .&(occursin.(Regex(a), dfa.CondVariables), occursin.(Regex(b), dfa.CondVariables)), dfa.CondStates .== cs12), :])
            dfz = vcat(dfz, dfa[.&(dfa.Target .== s, .&(occursin.(Regex(a), dfa.CondVariables), occursin.(Regex(c), dfa.CondVariables)), dfa.CondStates .== cs13), :])
            dfz = vcat(dfz, dfa[.&(dfa.Target .== s, .&(occursin.(Regex(a), dfa.CondVariables), occursin.(Regex(d), dfa.CondVariables)), dfa.CondStates .== cs14), :])
            dfz = vcat(dfz, dfa[.&(dfa.Target .== s, .&(occursin.(Regex(b), dfa.CondVariables), occursin.(Regex(c), dfa.CondVariables)), dfa.CondStates .== cs23), :])
            dfz = vcat(dfz, dfa[.&(dfa.Target .== s, .&(occursin.(Regex(b), dfa.CondVariables), occursin.(Regex(d), dfa.CondVariables)), dfa.CondStates .== cs24), :])
            dfz = vcat(dfz, dfa[.&(dfa.Target .== s, .&(occursin.(Regex(c), dfa.CondVariables), occursin.(Regex(d), dfa.CondVariables)), dfa.CondStates .== cs34), :])

            dfz = vcat(dfz, dfa[.&(dfa.Target .== s, .&(occursin.(Regex(a), dfa.CondVariables), occursin.(Regex(b), dfa.CondVariables), occursin.(Regex(c), dfa.CondVariables)), dfa.CondStates .== cs123), :])
            dfz = vcat(dfz, dfa[.&(dfa.Target .== s, .&(occursin.(Regex(a), dfa.CondVariables), occursin.(Regex(c), dfa.CondVariables), occursin.(Regex(d), dfa.CondVariables)), dfa.CondStates .== cs134), :])
            dfz = vcat(dfz, dfa[.&(dfa.Target .== s, .&(occursin.(Regex(b), dfa.CondVariables), occursin.(Regex(c), dfa.CondVariables), occursin.(Regex(d), dfa.CondVariables)), dfa.CondStates .== cs234), :])
            dfz = vcat(dfz, dfa[.&(dfa.Target .== s, .&(occursin.(Regex(a), dfa.CondVariables), occursin.(Regex(b), dfa.CondVariables), occursin.(Regex(c), dfa.CondVariables), occursin.(Regex(d), dfa.CondVariables)), dfa.CondStates .== cs1234), :])

        end

    end

    return dfz

end

