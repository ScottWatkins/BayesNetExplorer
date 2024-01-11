"""
    P(A=Astate,|B=Bstate,C=Cstate,...,N=n)

Parse a conditional probability expression.    
Return (targets, features, tt, ts, ff, fs) tuple.
Last four return values are inputs for bne.
"""
function queryparser(query::String)

    query = replace(query, " " => "")    # strip spaces

    if occursin(r"^P\(((.+)=(.+))\|((.+)=(.+))\)$", query)
        println("Query:$query")
        tf = split(query, r"\|")
        t = String.(split(tf[1][3:end], ","))
        f = String.(split(tf[2][1:end-1], ","))

        tt = [string.(split(i, "=")[1]) for i in t ]
        length(tt) == 1 ? tt=tt[1] : tt

        ts = [string.(split(i, "=")[2]) for i in t ]
        length(ts) == 1 ? ts=ts[1] : ts

        ff = [string.(split(i, "=")[1]) for i in f ]
        fs = [string.(split(i, "=")[2]) for i in f ]

        return(t, f, tt, ts, ff, fs)
    else
        error("\n\nBad conditional query: \"$query\"\nQuery format is \"P(A=Astate|B=Bstate,C=Cstate,...)\"\n")
    end

end
