# A stand alone script to parse the qbn condenser file.

names = []
OUT = open("dat.tmp", "w")
OUT2 = open("head.tmp", "w")
open(ARGS[1]) do f
    k = 0
    m = 0
 #   mflag = []
    
    for i in eachline(f)
        if occursin(r"^#A", i)

            b = split(i, "\t")
            c = replace(b[14], r"" => ",", count=length(b[14]))
            d = join([b[4],c], "")
#            TRU = parse(Int, b[8])
#            d = join([d,TRU], ",")

#            if TRU == 2
#                m = m + 1
#            end        

            k = k + 1
            id = "ID" * string(k) * ","
            println(OUT, id, d)
        end

        if occursin(r"^#O", i)
            o = split(i, r"\s+")
            push!(names, "ID", replace(o[3], "," => "_"))
        end
        
        if occursin(r"^#B", i)
            x = split(i, "\t")
            push!(names, rstrip(x[4]))
        end

    end

    nm = join(names, ",")
#    nm = nm * ",mflag"
    println(OUT2, nm)

#    println("Found $m individuals with missing indicator for condenser file (TRU=2)")
#    println("This condenser file specific missingness is now indicated in the mflag column")

end

close(OUT)
close(OUT2)

run(pipeline(`cat head.tmp dat.tmp`, stdout="qbn_matrix.csv"))
run(`rm head.tmp dat.tmp`)
println("Output written to qbn_matrix.csv")
