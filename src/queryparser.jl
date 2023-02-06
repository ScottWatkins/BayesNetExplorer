"""
    P(A=a,|B=b,C=c,...,N=n)

Parse a conditional probability expression.    
Return (target, feature) tuple.
"""
function queryparser(query)
    if occursin(r"P(.+)", query)
        tf = split(query, r"\|")
        t = String.(split(tf[1][3:end], ","))
        f = String.(split(tf[2][1:end-1], ","))
        return(t, f)
    else
        error("\n\nInput query must be specified as P(X=xi|Y=yi) not $query \n\n")
    end
    
end
