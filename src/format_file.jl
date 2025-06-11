"""
    df, ids, features = format_file(datafile; datacols=[], delim=",", minfreq=0.0, recode_bool=false)

Create and format variables to be used in the network from a user provided data table. Write data and header files to disk (BN.data, BN.header). Write a variable map file to disk (recode.map). The user input table *must* be text based unless it contains strictly boolean variables coded as 0/1 or true/false.

The input data must be a full matrix with observations in rows and features in columns. Observations and features must be fully labeled.  Features must be discrete or quantized variables. 


Options: 

    datacols        columns to include         [:IDS, :VAR1, :VAR2]  
    recode_bool     input is strictly boolean  [false]  
    delim           set delimiter              [","]  
    minfreq         min variable state freq    [0.0]

Example text dataset: 

ID,SIZE,FLIGHT,HABITAT\\
X1,big,N,SAVANA \\
X2,med,N,JUNGLE \\
S4,sml,Y,FOREST \\
D7,med,Y,RIPARIAN \\


Example Boolean dataset:

ID,V1,V2,V3 \\
ID1,1,0,true \\
ID4,0,0,false \\
ID8,1,0,false \\
ID9,1,1,true \\

"""
function format_file(infile::Union{String,DataFrame}; datacols::Array=[], delim::Union{String,Char}=",", minfreq::Float64=0.0, recode_bool::Bool=false, recode12::Bool=false, datafile="BN.data", headerfile="BN.header")

    if typeof(infile) == String
        if !isfile(infile)
            error("\n\nCould not find infile: $infile\n\n")
        end
    end

    df, rmap = recoder(infile; recode_bool=recode_bool, recode12=recode12, delim=delim)
println(names(df))
    ids = df[:,1]
    dfc = size(df,2)

    BN = df
    col1 = Symbol(names(df)[1])
    #println("==>", col1)
    if length(datacols) > 0

        if in(1, datacols) || in(Symbol(names(df)[1]), datacols)
            println("Selected data includes the labels column...")
        else
            error("\n\nRemember to always include the first column, don't exceed column the total columns in this dataframe ($dfc), and use recode_bool=true if the input file is strictly 0/1 variables.\n\n")
        end
        
        BN = select(df, datacols)

    end

    if minfreq > 0.0

        dfcleaned, df_freq, labels = clean_dfvars_by_frequency(BN, minfreq, header=true, nolabels=false)
        BN = insertcols!(dfcleaned, 1, col1 => ids)

    end            
    
    if recode_bool == false

        BN = BN[!, 2:end]
        CSV.write(datafile, BN, delim=" ", header=false)
        dorc = [names(BN), length.(unique.(eachcol(BN))), repeat(["D"], inner=size(BN,2) )]

    elseif  recode_bool == true

        BN = BN[!, 2:end]
        BN[:,:] = replace.(BN[:,:], "true" => "2", "false" => "1")
        CSV.write(datafile, BN, delim=" ", header=false)
        dorc = [names(BN), length.(unique.(eachcol(BN))), repeat(["D"], inner=size(BN,2) )]

    else
        
        cols = names(df)[2:end]
        BN = df[:, 2:end]
        CSV.write(datafile, BN, delim=" ", header=false)
        dorc = [names(BN), length.(unique.(eachcol(BN))), repeat(["D"], inner=size(BN,2) )]

    end
    
    OUT = open(headerfile, "w+")

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
        rm(headerfile)
        rm(datafile)
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
