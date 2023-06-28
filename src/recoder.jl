"""
    r, rmap = recoder(d::Union{String,DataFrame}; recode_bool=false, delim=",")

Recode text data from a file or dataframe to a numerical or boolean dataframe. Recoded values are returned in the rmap dictionary. Input data must be a full matrix with row and column labels.

Options:
    recode_bool    Convert 0/1 dataframe to true/false and write to disk
    recode12       Convert 0/1 dataframe to 1/2 dataframe

"""
function recoder(d::Union{String,DataFrame}; recode12=false, recode_bool=false, delim=delim)

    if typeof(d) == String && isfile(d) == true
        println("Reading data from file: $d")
        d = CSV.read(d, DataFrame, delim=delim, types=String)
    end
    
    col1name = Symbol(names(d)[1])

    if recode_bool == true

        db = Bool.(d[!, 2:end])
        d = hcat(d[:,1], db)
        rename!(d, :x1 => col1name)
        CSV.write("recoded.bool.csv", d, delim=",")
        println("Wrote file recoded.bool.csv to disk.")
        d = CSV.read("recoded.bool.csv", DataFrame, delim=",", types=String)
        
    elseif  recode12 == true

        nl = names(d)[1]
        l = d[:,1]
        d = d[:, 2:end] .+ 1
        
        insertcols!(d, 1, "$nl" => l)
        CSV.write("d.tmp", d, delim=",")
        d = CSV.read("d.tmp", DataFrame, delim=",", types=String)
        rm("d.tmp")

    else
        
        if in(0, unique(d[:,2]) ) 
            if !in(2, unique(d[:,2]) )
                println("Found values: ", unique(d[:,2]), " for first variable.")
                error("\n\nInput data appears to be Boolean. Use recode_bool=true to format a strictly 0/1 or true/false data set. Alternatively, boolean column data can be reformated to  strings (e.g. Y/N). Except for ids (column 1), input data must be uniformly Boolean or String type. \n\n")
            end

            if maximum(d[!, 2:5]) > 1
                error("Please code data with multistate variables as strings, e.g. big/med/small ")
            end
            
        end
    end

    printstyled("Recoder matrix input: $(size(d)).\nRecoding and mapping all variables...\n", color=:blue)
    
    rmap = Dict()
    vcount = []

    for c in 2:size(d, 2)
        
        d[!,c] = string.(d[!,c]) #stringify all
        
        u = sort(unique(d[!, c])) #false is 1, true is 2
        push!(vcount, length(u))

        k = 1
        
        for i in 1:length(u)
            
            if recode12 == false && recode_bool == false
                d[:,c] = replace!(d[:,c], u[i] => string.(k))
            end
            
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
