"""
    P(A=a,|B=b,C=c,...,N=n)

Parse a conditional probability expression.    
Return (target, feature) tuple.
"""
function queryparser(query)

    query = replace(query, " " => "")    # strip spaces

    if occursin(r"^P\(((.+)=(.+))\|((.+)=(.+))\)$", query)
        println("Query:$query")
        tf = split(query, r"\|")
        t = String.(split(tf[1][3:end], ","))
        f = String.(split(tf[2][1:end-1], ","))
        return(t, f)
    else
        error("\n\nBad conditional query: \"$query\"\nQuery format is \"P(A=Astate|B=Bstate,C=Cstate,...)\"\n")
    end

end
