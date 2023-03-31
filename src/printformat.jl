"""
    printformat(n; fmt=f)

Format a variable with printf. The printf format
string is passed as a variable.
"""
function printformat(n; fmt="")
    Printf.format(Printf.Format(fmt), n)
end
