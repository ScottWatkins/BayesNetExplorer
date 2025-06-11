"""
    phenomaker(df; newpheno="", options...)

Create a new phenotype by combining variables using AND/OR conditions and add the new variable to the dataframe.

Input dataframe column 1 MUST contain labels such as IDs or sample names.

    Required
        newpheno        name for new composite phenotype       ""

    Options
        mergeall        merge expression or phenos in a file   ""
                        
        expression      \"\"\"AND/OR expression\"\"\"                \"\"\" \"\"\"

        normalizenames  normalize all input names              true
                        (false not recommended)

For expressions, specify combinations of variables in a *triple quoted* expression with varible states in single quote ("0" or "1"). Data in the input dataframe should be coded as 1 = presense of a phenotype and 0 = absense of a phenotype. Dataframe column types must be strings. Columns names will be uppercased and special characters removed unless normalizenames is set to false.

In the expression string, .& means AND while .| means OR. Combinations must nest into a single expression in the form .&(...) or .|(...). The expression should not duplicate any column names.

Examples
--------

1. Merge T2D and BMI40 phenotyes into a new composite phenotype

   phenomaker(df, newpheno=T2DhighBMI, mergeall=\"\"\"T2D,BMI40\"\"\" )

2. Merge all phenotype using a list of phenotypes in a file.

   phenomaker(df, newpheno=T2DhighBMI, mergeall="phenolistfile.txt" )

3. Merge T2D AND BMI40 OR high glucose using an expression

   phenomaker(df, newpheno=T2DBMIallGLU, expression=\"\"\".|(GLU="1", .&(T2D="1",BMI40="1"))\"\"\" )


Additional expression examples (white space is optional)
--------------------------------------------------------

\"\"\".&(pheno1 = "1", pheno2 = "1")\"\"\"     intersect

\"\"\".|(pheno1 = "1", pheno2 = "1")\"\"\"     union

\"\"\".|( .&(pheno1 = "1", pheno2 = "0"), (pheno3 = "1") )\"\"\"      intersect 1,2 then union 3

\"\"\".|( .&(pheno1="1", pheno2="0"), .|(pheno3="1", pheno4="1") )\"\"\"  1 and 2 OR 3 and 4

\"\"\".|(.&(pheno1 = "1", pheno2 = "0")), .&(pheno3 = "1")\"\"\"   WRONG, non-nesting!

"""
function phenomaker(df::DataFrame; expression::String="", newpheno::Union{String,Symbol}="", mergeall::String="", normalizenames=true)

    ucn = uppercase.(names(df))
    df = rename!(df, uppercase.(names(df)))
    
    if length(mergeall) > 0

        if length(expression) > 0
            error("\n\nPlease use either the mergeall or expression keyword. Use mergeall only when all phenotypes are to be combined into a single composite phenotype.\n\n")
        end
        
        df2 =  targetmaker(df, mergeall, newpheno )

        return df2

    else

        global dfg = df

        if occursin.(r".==", expression)
            error(".== should be = only.")
        elseif occursin(r"^.[&|\|]", expression)
            println("Processing $expression ...")
        else
            error("\n\nInput expression should begin with .& or .|\n
    Example:  \"\"\".&(pheno1=\"1\", pheno2=\"0\")\"\"\"
    Example:  \"\"\".|(pheno1=\"0\", .&(pheno2=\"1\", pheno3=\"1\"))\"\"\" \n\n")
        end

        a = expression
        e = expression
    
        if normalizenames == true
            a = replace.(a, r".[&|\|]\(" => "")
            a = replace.(a, ")" => "")
            a = uppercase.(replace.(a, r"[.|\$/@^]"  => "_"))
        end

        b = string.(strip.(split(a, r"[(=,]")))
        b = b[occursin.(r"^[A-Z]+", b)]  # clean user col names
        
        c = findall(in(b), names(df))    # column idx for targets

        if length(b) == length(c)
            println("All input column names matched...")
        else
            cn = names(df)
            for i in b
                ir = Regex("^" * i * "\$")
                if sum(occursin.(ir, cn)) > 0
                    println("Column $i ... found.")
                else
                    printstyled("Column $i ... not found.\n", color=:yellow)
                end
            end
            error("\n\nNot all columns found or duplicate columns requested...\nPlease check list above. Variable names are case sensitive.
    Option normalizenames was set to $normalizenames\n\n")
        end
        
        if length(collect(eachmatch(r"\(", expression))) == length(collect(eachmatch(r"\)", expression)))
            println("Expression check ... OK")
        else
            error("\n\nMismatched parentheses in expression!\n\n")
        end

        expression = replace(expression, r"\s+" => "")
        expression = replace(expression, "=" => ".==") 
        expression = replace(expression, r"(,)([A-Za-z])" => s"\1dfg.\2" )
        expression = replace(expression, r".\&\(([A-Za-z])" => s".&(dfg.\1")
        expression = replace(expression, r".\|\(([A-Za-z])" => s".|(dfg.\1")

        q = "dfg[" * expression * " , 1]"
        q = Meta.parse(q)
        eq =eval(q)    # array of ids matching user expression
        
        c = string.(Int.(zeros(size(df,1))))
        idxm = []

        for i in eq
            idxm = findall(in(eq), df[:,1])
            c[idxm] .= "1"
        end

        insertcols!(df, 2, newpheno => c)
        cs = countmap(df[!,2])
        println("Created new variable $newpheno ...")
        println("State\tCount")
        println("-"^14)

        for (k,v) in cs
            println(k, "\t", v)
        end

        println("-"^14)
    
        return df
    
    end

end
