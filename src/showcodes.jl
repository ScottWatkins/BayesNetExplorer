"""
    showcodes()

Display the data variables and their coded values.

"""
function showcodes()
    dfc = CSV.read("recoded.map", DataFrame, types=String)
    dfc = rename(dfc, [:VAR, :VAR_STATE, :VAR_CODE])
    println(dfc)
end
