"""

    ivarsall_set = create_powerset(df_freq::DataFrame)

Create a powerset of combinations for two arrays. Read arrays from a dataframe. Combine arrays and create combinations.
"""
function create_powerset(df_freq::DataFrame)

    #Create merged variables and states (A, B, 1, 2 --> A1, B1, A2, B2)                                                                         
    ivarsall = [df_freq.feature[i] * ":" * string(df_freq.value[i]) for i in eachindex(df_freq.feature)]
    
    #Selector to remove the target from all states                                                                                              
    rmf = [!i for i in occursin.(f, ivarsall)]
    ivarsall = ivarsall[ rmf ]
    
    #Create all combinations of the merged data                                                                                                 
    ivarsall_set = collect(powerset(ivarsall))

    return ivarsall_set
end

