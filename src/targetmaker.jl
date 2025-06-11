"""
    targetmaker(df, targets, newcolname; normalizenames=true)

Merge specified 0/1 string variables in a dataframe into a single new variable, and add the new variable to column two in the dataframe.

The data type for every column should be String. Values other than "1" or "0" are ignored.

Specify selected target columns in a string  "target1,target2, ..." or in a text file, one column name per line. Columns names will automatically be uppercased and special characters removed unless normalizenames is set to false.

Examples
--------
targetmaker(df, "ICD10DX.J45.21,ICD10DX.J45.909", "Asthma")

targetmaker(df, "asthma.codes.txt", "Asthma")

"""
function targetmaker(df::DataFrame, targets::String, cname::String; normalizenames=true)

    if (length(targets) < 100) && isfile(targets)
        a = strip.(split(read(targets, String), "\n"))
        null = pop!(a)
    else
        a = strip.(string.(split(targets, ",")))
    end

    if normalizenames == true
        a = uppercase.(replace.(a, r"[.|\$/@^)(]"  => "_"))
    end
    
    b = findall(in(a), names(df))    # column idx for targets
    
    if length(a) == length(b)
        println("All input column names matched...")
    else
        cn = names(df)
        for i in a
            ir = Regex("^" * i * "\$")
            if sum(occursin.(ir, cn)) > 0
                println("Column $i ... found.")
            else
                printstyled("Column $i ... not found.\n", color=:yellow)
            end
        end
        error("\n\nNot all columns found...\nPlease check list above.\n\n")
    end
 
    c = Int.(zeros(size(df,1),1))    # new array

    for i in b
        idxm = findall(df[:,i] .== "1")
        c[idxm] .= 1
    end

    c = vec(string.(c))

    cs = countmap(c)

    println("Created new variable $cname ...")
    println("State\tCount")
    println("-"^14)
    for (k,v) in cs
        println(k, "\t", v)
    end
    println("-"^14)
    
    insertcols!(df, 2, cname => c)

    return df
    
end

