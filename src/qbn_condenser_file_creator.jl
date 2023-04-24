"""

    qbncfc(df, target)

Create a qbn infile from a dataframe.
Input dataframe must be 0/1 where 1 represents a test condition/disease.
"""
function qbncfc(df::DataFrame, target::String)

    OUT = open("condenser.out", "w")
    
    I = rpad("I",7); PID=rpad("PID", 64)
    pid = randstring(['a':'z'; '0':'9'], 64)
    pid = bytes2hex(sha256(pid))
    f = (sum( df[!, Symbol("$target")] ) ) / size(df,1)

    println(OUT, "#O $f    O($target)")
    println(OUT, "##\tI\t$PID\tO\tBEG\tEND\tDUR\tTRU\tAGE\tVST\tLST\tEXP\tTAG\tIO_STR" )

    for i in 1:size(df,1)
        pid = randstring(['a':'z'; '0':'9'], 64)
        id = i - 1
        t = df[i, Symbol("$target")]
        
        BEG, END, DUR = 999, 999, 999;
        t == 1 ? TRU = 1 : TRU = 0
        AGE = rand([20:85...] ,1)[1]
        VST = rand([2:10...], 1)[1]
        LST = rand([1995:2021...], 1)[1]
        EXP = VST*AGE
#        EXP = rand([1000:2000...], 1)[1]
        t == 1 ? TAG = "0(1,0)" : TAG = "0(0,0)"
        IO_STR = join(df[i, 3:end], "")
        
        println(OUT, "#A\t", id, "\t", pid, "\t", t,  "\t", BEG, "\t", END, "\t", DUR, "\t", TRU, "\t", AGE, "\t", VST, "\t", LST, "\t", EXP, "\t", TAG, "\t", IO_STR)
    end
    println(OUT, "##\tI\t$PID\tO\tBEG\tEND\tDUR\tTRU\tAGE\tVST\tLST\tEXP\tTAG\tIO_STR" )

    println(OUT, "##      OFFSET  FRQ             TID/BIN                         DESCRIPTION")

    for i in 3:size(df,2)
        col = names(df)[i]
        offset = i - 3
        frq = rpad(round(sum(df[!,Symbol("$col")]) / size(df,1), digits=5), 10)
        println(OUT, "#B\t", offset, "\t", frq, "\t", rpad(names(df)[i], 25), "\t", "Trait $offset")
    end
    
    println(OUT, "##      OFFSET  FRQ             TID/BIN                         DESCRIPTION")
    println(OUT, "## NUMBER OF PROBANDS:",size(df,1))
    println(OUT, "##This condenser output file was generated from a datatable.")
    println(OUT, "##2020-01-01 01:01:01.000001")
    println(OUT, "##EOF")


    close(OUT)

    printstyled("WARNING: the output is currently not working for qbn input (Feb 2023)", color=:red)

end
    
