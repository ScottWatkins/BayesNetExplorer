"""
    merge_bne_cpq(bnefile, cpqfile; addheader=true, outfile="merged.csv")

A utility function to merge identical queries from network bne runs and simple conditional probabilities from cpq runs. The queries should be the same. If not, the cpqfile queries are dropped. File will be headered with appropriate header files.
"""
function merge_bne_cpq(bnefile::String, cpqfile::String; addheader=true, outfile="merged.csv")

    if isfile(bnefile) && isfile(cpqfile)

        if addheader == true
            run(`bash -c "cat bootheader.txt $bnefile > b; mv b $bnefile"`)
            run(`bash -c "cat header.txt $cpqfile > c; mv c $cpqfile"`)
            println("Added headers to files...")
        end
        
        dfb = CSV.read(bnefile, DataFrame, delim="\t", comment="#")
        dfc = CSV.read(cpqfile, DataFrame, delim="\t", comment="#")

        dfz = leftjoin(dfb, dfc, on=:Query)
        println(dfz)
        
        if isequal(missing, dfz.Binomial_pval)
            println("A complete list of all binomial p-values is needed for the BH correction.")
        else
            println("Found Binomial_pvals...")
            bjc = benjhoc(dfz.Binomial_pval)
            insertcols!(dfz, 26, :Binomial_pval_BH_corrected => bjc)
            println("Added Benjamini-Hochberg corrected binomial values...")
        end

        println("Merged results written to disk file: mergedresults.csv")

        CSV.write(outfile, dfz)
        
        return dfz

    end        
end
