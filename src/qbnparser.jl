# A stand alone script to parse the qbn condenser file.

names = []
OUT = open("dat.tmp", "w")
OUT2 = open("head.tmp", "w")
open(ARGS[1]) do f
    k = 0
    for i in eachline(f)
        if occursin(r"^#A", i)
            b = split(i, "\t")
            c = replace(b[14], r"" => ",", count=length(b[14]))
            d = join([b[4],c], "")
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
    println(OUT2, nm)
end

close(OUT)
close(OUT2)

run(pipeline(`cat head.tmp dat.tmp`, stdout="qbn_matrix.csv"))
run(`rm head.tmp dat.tmp`)
println("Output written to qbn_matrix.csv")
