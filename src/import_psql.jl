"""

    df = import_psql(dbname, user)

Connect and retrieve data from PostgreSQL.

    Options:
        password    user password if needed    ""
        host        hostname                   "localhost"
        port        connection port            5432
        
Dependencies: LibPG, Tables, DataFrames, and a running PostgreSQL instance.
"""
function import_psql(dbname::String, user::String; password::String="", host="localhost", port::Int64=5432)

    
    conn = LibPQ.Connection("dbname=$dbname user=$user password=$password host=$host port=$port")

    println("Tables available in $dbname for query are:")
    tlist = "SELECT schemaname, tablename FROM pg_catalog.pg_tables WHERE tablename NOT LIKE 'pg_%' AND tablename NOT LIKE 'sql_%';"
    dft = execute(conn, tlist) |> DataFrame
    println(dft)

    println("List columns for a table (y/n)?")
    q = chomp(readline())
        
    while occursin(r"Y|y", q)
        println("Enter table name:")
        tn = readline()
        cq = "SELECT table_name , column_name FROM information_schema.columns WHERE table_name = '$tn';"
        dfcq = execute(conn, cq) |> DataFrame
        println(dfcq)
        println("List columns for another table (y/n)?")
        q = chomp(readline())
    end

    println("Type your SQL statment then <enter>")
    sql_in = readline()
    df = execute(conn, sql_in) |> DataFrame

    close(conn)
    
    return df
    
end
