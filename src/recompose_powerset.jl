"""
g, gs = recompose_powerset(psets)

Take individual elements of a powerset array and decompose
to two separate arrays.

If you need a powerset that is a combination of two arrays of elements, ([A,B,C] and [1,2]), first combine the elements indat = [A1, A2, B1, B2, C1, C2] and run x = collect(powerset(indat)). x will be [[], [A1], [A2] ... [A1, A2, C2] etc. This function then decomposes any x[i] to a pair of arrays g, gs as in g=[A, A, C] and gs = [1, 2, 2].
"""
function recompose_powerset(psets)
    g = []
    gs = []    
    for i in eachindex(psets)
        decomp = split.(psets[i], ":")
        push!(g, decomp[1])
        push!(gs, decomp[2])
    end
    #print(g)
    #print(gs)
    return(g,gs)
end

