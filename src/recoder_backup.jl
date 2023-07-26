"""A
    r, rmap = recoder(d::Union{String,DataFrame}; boolmat=false, delim=",")

Recode text data from a file or dataframe to a numerical or boolean dataframe. Recoded values are returned in the rmap dictionary. Input data must be a full matrix with row and column labels.

Options:
    boolmat        Convert 0/1 dataframe to true/false and write to disk
    recode12       Convert 0/1 dataframe to 1/2 dataframe

"""
function recoder(d::Union{String,DataFrame}; recode12=false, boolmat=false, delim=",")

    if typeof(d) == String && isfile(d) == true
        println("Reading data from file: $d")
        d = CSV.read(d, DataFrame, delim=delim)
    end

    if boolmat == true
        db = Bool.(d[!, 2:end])
        d = hcat(d[!,1], db)
        rename!(d, :x1 => "IID")
        CSV.write("recoded.bool.csv", d, delim=",")
        println("Wrote file recoded.bool.csv to disk.")
        d = CSV.read("recoded.bool.csv", DataFrame, delim=",", types=String)
    end
    
    if recode12 == true
        nl = names(d)[1]
        l = d[!,1]
        d = d[!, 2:end] .+ 1
        insertcols!(d, 1, "$nl" => l)
        CSV.write("d.tmp", d, delim=",")
        d = CSV.read("d.tmp", DataFrame, delim=",", types=String)
        rm("d.tmp")
    end

    #println(d)

    rmap = Dict()
    vcount = []

    for c in 2:size(d, 2)

        u = unique(d[!, c])
        push!(vcount, length(u))

        k = 1

        for i in 1:length(u)
            d[:,c] = replace!(d[:,c], u[i] => string.(k))
            m = join([ names(d)[c], u[i] ], ":")
            rmap[m] = string(k)
            k += 1
        end

    end


    CSV.write("recoded.out", d, delim=",")
    OUT = open("recoded.map", "w")
    println(OUT, "feature,state,numstate")
    for i in sort(collect(keys(rmap)))
        (k1,k2) = split(i, ":") 
        println(OUT, k1, ",", k2, ",", rmap[i])     
    end
    close(OUT)
    println("Wrote recode data and recode map to recoded.out and recoded.map.")

    return d, rmap

end
