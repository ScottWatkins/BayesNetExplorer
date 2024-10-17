"""
    P(A=Astate,|B=Bstate,C=Cstate,...,N=n)

Parse a conditional probability expression.    
Return (targets, features, tt, ts, ff, fs) tuple.

"""
function queryparser(query::String)

    query = replace(query, " " => "")    # strip spaces

    if occursin(r"^P\((.+)=(.+)\|(.+)=(.+)\)$", query) #basic check

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

        error("\n\nBad query: \"$query\"\n\nPlease use the standard query format, for example:\n \"P(A=a|B=b,C=c, ...)\" or \"P(A=a,B=b|C=c,D=d, ...)\" \n and limit targets to no more than three variables\n")

    end

end
