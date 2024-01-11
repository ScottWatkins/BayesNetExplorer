"""
    mergedat(outfilename)

Combine all permutation runs from each target.
"""
function mergedat(filename)
    run(`bash -c "head -n1 RRtable_MORT.csv > header;
    cat RRtable_*.csv >RRtable_ALL.csv;
    sed '/^Target/d' RRtable_ALL.csv >j;
    cat header j >k;
    mv k RRtable_ALL.csv;
    rm -f j;
    cp RRtable_ALL.csv $filename;" `)
end
