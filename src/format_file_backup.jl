"""
    makeBNstructfiles(datain; datacols=[], delim=",", clean=false, minfreq=0.0, recode_bool=false, recode12=false)

Write a bnstruct input data matrix and header file from a data source. 

Input data can be a disk file or dataframe and must be a full matrix of observations in rows and features in columns. Observations and features must be fully labelled.  Features must be discrete or quantized variables. 

Options:

datacols       specify a subset of data columns with vector array [].
recode_bool    recode 0/1 input to boolean and 1/2 variables [false]
delim          set delimiter for data file input [","]
clean          clean the data [false]
minfreq        min frequency for any variable state [0.0]. Use with clean.
recode12       recode the output as 1 and 2; rarely used

"""
function makeBNstructfiles(infile::Union{String,DataFrame}; datacols::Array=[], delim::Union{String,Char}=",", clean::Bool=false, minfreq::Float64=0.0, recode_bool::Bool=false, recode12::Bool=false)

    df, rmap = recoder(infile, recode_bool=recode_bool, recode12=recode12, delim=delim)

    if length(datacols) > 0 && clean == true
        error("Manual column selection is not available with frequency-based variable selection.")
    end

    if clean == false && minfreq > 0.0
        error("Please set clean to true to enable data filtering options.")
    end    
    
    if clean == true
        
        if minfreq > 0.0
            obscol = df[:,1]
            dfcleaned, df_freq, labels = clean_dfvars_by_frequency(df, minfreq, header=true, nolabels=false);
            df = insertcols!(dfcleaned, 1, :IID => obscol)
        end        
    end
    
    if recode_bool == false
        
        if length(datacols) > 0
            
            if in(1, datacols)
                error("\n\nThe datacols vector should not include the id column 1.\n")
            end

            if maximum(datacols) > size(df,2)
                error("\n\nMaximum value in datacols vector exceeds number of columns in data.\n")
            end
            
            pushfirst!(datacols, 1)
            BN = df[!, datacols]

        else    
            BN = df
        end
        
        #println(BN)
        BN = BN[!, 2:end]
        CSV.write("BN.data", BN, delim=" ", header=false)
        dorc = [names(BN), length.(unique.(eachcol(BN))), repeat(["D"], inner=size(BN,2) )]

    elseif  recode_bool == true

        if length(datacols) > 0
            BN = select(df, datacols)
        else
            BN = df
        end
        
        BN[:,:] = replace.(BN[:,:], "true" => "2", "false" => "1")
        println(BN)
        BN = BN[!,2:end]
        CSV.write("BN.data", BN, delim=" ", header=false)
        dorc = [names(BN), length.(unique.(eachcol(BN))), repeat(["D"], inner=size(BN,2) )]

    else
        
        cols = names(df)[2:end]
        BN = df[:, 2:end]
        CSV.write("BN.data", BN, delim=" ", header=false)
        dorc = [names(BN), length.(unique.(eachcol(BN))), repeat(["D"], inner=size(BN,2) )]

    end
    
    OUT = open("BN.header", "w+")

    qd = join(["\"" *  i *  "\""  for i in names(BN)], " ")
    println(OUT, qd)

    c = join( [ i for i in length.(unique.(eachcol(BN))) ], " ")
    println(OUT, c)

    s = join( [ i for i in dorc[3] ], " ")

    println(OUT, s)

    close(OUT)

    if size(BN) !== size(dropmissing(BN))

        md = 1
        miss = size(BN, 1) - size(dropmissing(BN), 1)
        printstyled("INFO: Detected $miss rows of missing data. Please use the impute_dataframe.jl function to impute the missing data!\n", color=:yellow)
        rm("BN.header")
        rm("BN.data")
        error("Input data has missing data.")

    else

        md = 0
        printstyled("INFO: Full matrix of data detected $(size(BN)).\n", color=:green)
        printstyled("INFO: Using columns: $datacols\n", color=:green)
        
    end

    println("Wrote BN.data and BN.header files to disk.")

    return BN
end
