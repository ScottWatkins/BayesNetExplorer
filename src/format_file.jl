"""
    df, ids, features = format_file(userdata; datacols=[], delim=",", clean=false, minfreq=0.0, recode_bool=false, recode12=false)

Format input data and header files from a user provided data table. Write data and header files to disk (BN.data, BN.header). Write variable map file to disk (recode.map). 

Input data must be a full matrix with observations in rows and features in columns. Observations and features must be fully labeled.  Features must be discrete or quantized variables. The user input table must be text based unless it contains strictly boolean columns coded as 0/1 or true/false.

Input examples

ID,SIZE,FLIGHT,HABITAT \\
X1,big,N,SAVANA \\
X2,med,N,JUNGLE \\
S4,sml,Y,FOREST \\
D7,med,Y,RIPARIAN \\


ID,V1,V2,V3 \\
1,1,0,true \\
4,0,0,false \\
8,1,0,false \\
9,1,1,true \\


Options:

datacols: specify a subset of data columns with vector array []  \\
recode_bool: recode 0/1 input to boolean and 1/2 variables [false]  \\
delim: set delimiter for data file input [","]  \\
minfreqmin: frequency for any variable state [0.0] Use with clean  \\

"""
function format_file(infile::Union{String,DataFrame}; datacols::Array=[], delim::Union{String,Char}=",", minfreq::Float64=0.0, recode_bool::Bool=false, recode12::Bool=false)

    df, rmap = recoder(infile, recode_bool=recode_bool, recode12=recode12, delim=delim)

    ids = df[:,1]
    dfc = size(df,2)

    BN = df

    
    if length(datacols) > 0

        if !in(1, datacols) 
            error("Please always include the labels (column 1) and don't exceed column $dfc !\n\n")
        end
        
        if maximum(datacols) > dfc
            error("Data column values should not exceed $dfc")
        end
        
        
        BN = select(df, datacols)

    end

    if minfreq > 0.0

        dfcleaned, df_freq, labels = clean_dfvars_by_frequency(BN, minfreq, header=true, nolabels=false)
        BN = insertcols!(dfcleaned, 1, :IID => ids)

    end            
    
    if recode_bool == false

        BN = BN[!, 2:end]
        CSV.write("BN.data", BN, delim=" ", header=false)
        dorc = [names(BN), length.(unique.(eachcol(BN))), repeat(["D"], inner=size(BN,2) )]

    elseif  recode_bool == true

        BN = BN[!, 2:end]
        BN[:,:] = replace.(BN[:,:], "true" => "2", "false" => "1")
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
        printstyled("INFO: Detected $miss rows with missing data. Please use the impute_dataframe function to impute the missing data!\n", color=:yellow)
        rm("BN.header")
        rm("BN.data")
        error("Input data has missing data.")

    else

        md = 0
        printstyled("INFO: Full matrix of data detected $(size(BN)).\n", color=:green)
        printstyled("INFO: Using columns: $datacols\n", color=:green)
        
    end

    println("Wrote BN.data and BN.header files to disk.")

    ids = df[!, 1]
    features = names(BN)
    
    return BN, ids, features

end
