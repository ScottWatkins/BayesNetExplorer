"""
    querywriter(RRtable)

Returns an array of probability queries from a RRtable for input into bnemle.
 
"""
function querywriter(rrtable::DataFrame)

    qa = []

    for i in eachrow(rrtable)
           cv = split(i[9], ",")
           cs = split(i[11], ",")
           condvars=""

        for j in 1:length(cv)
            condvars = condvars * cv[j] * "=" * cs[j] * ","
        end

        condvars = condvars[1:end-1]
        push!(qa, "P(" * i[1]* "=2|" * condvars * ")")

    end

    return qa

end
