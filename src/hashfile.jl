"""
    hashfile("inputfile", keycol, [valcols])

Create a dictionary from a tsv file with some column of the file as the lookup key and values from columns specified in an array. Comment lines starting with # are always omitted. For multiple keys, values are pushed to one value array.

kwargs: delim="\t", omit="", limit=0

    delim      set the file delimiter
    omit       omit the line if it contains this string
    omit2      a optional second exlusion string
    limit      max fields after split (e.g. 2 -> split on first delim only)
"""
function hashfile(file::String, keycolumn::Int, valcols::Array; delim="\t", omit="^\$", omit2="a^", limit=0)

    fhash = Dict{String,Array}()
    d = Regex(delim)
    o = Regex(omit)
    o2 = Regex(omit2)
    h = Regex("^#")
    
    open(file) do f

        for i in eachline(f)

            if occursin(h, i) || occursin(o, i) || occursin(o2, i)
                continue
            end
            
            line = split(i, d, limit=limit)

            if haskey(fhash, line[keycolumn])
                    push!(fhash[line[keycolumn]],  join(line[valcols], delim) )
            else
                fhash[line[keycolumn]] = [ join(line[valcols], delim) ] 
            end

        end
    end
    
    rhash = Dict{String,String}()

    for (k,v) in fhash
        sv = join(unique(v), ";")
        k = uppercase.(replace.(k, r"[.|\$/@^)(]"  => "_"))
        rhash[k] = sv
    end
    
    return rhash

end
