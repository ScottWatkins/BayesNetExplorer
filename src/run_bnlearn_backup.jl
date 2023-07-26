"""
    netout, plotout = run_bnlearn("datafile"; algo="hc", method="bayes", score="bic")

    Construct a Bayesian network using bnlearn. Input is a full csv dataset with labels, header, and no missing data. This function is a wrapper for R's bnlearn. R must be installed.

    Network learning methods:\\
    hc (greedy hill climbing method)\\
    tabu (greedy tabu search)\\
    mmhc (max-min hill-climbing heuristic)\\
    h2pc (hybrid HPC)\\
\\
    Parameter methods:\\
    bayes (default)\\
    mle (maximum likelihood estimator)
\\
    Network scoring methods:\\
    aic (Akaike Information Criterion)\\
    bic (Bayesian Information Criterion)\\
"""
function run_bnlearn(data::Union{DataFrame, String}; algo="hc", method="bayes", score="bic")
    
    if typeof(data) == String && isfile(data)
        df = CSV.read(data, DataFrame, delim=",", header=true)
        df = Bool.(df[:, 2:end])
        CSV.write("BNinfile.csv", df, delim=",", header=true)
    elseif typeof(data) == DataFrame
        df = data
        CSV.write("BNinfile.csv", df, delim=",", header=true)
        df = Bool.(df[:, 2:end])
        CSV.write("BNinfile.csv", df, delim=",", header=true)
        data = "BNinfile.csv"
    else
        error("A problem with the input file has occurred. Please check the input.")
    end    
    
    df_size = size(df)
    header = names(df)
    
    if df_size !== size(dropmissing(df))
        error("Input data set has missing data. Missing data must be imputed.")
    else
        print("Loaded data matrix.")
        println(df_size)
    end
    
    print("Loading R libraries...")
    R"library(bnlearn)"; R"library(Rgraphviz)";
    R"library(igraph)"; R"library(dplyr)"; R"library(stringr)";
    R"library(visNetwork)"; R"library(networkD3)"
    println("done.")
    
    println("Input data file: $data")
    println("Input header file: $header")

    print("Reading disk file...")
    R"dataset <- read.csv('BNinfile.csv', header = T, stringsAsFactors=TRUE)"
    println("done.")


    print("Learning net structure ($algo)...")
    if algo == "hc"
        R"dagout = hc(dataset)"    ##R obj, dag 
    elseif algo == "tabu"
        R"dagout = tabu(dataset)"  ##R obj, dag 
    elseif algo == "mmhc"
        R"dagout = mmhc(dataset)"  ##R obj, dag 
    elseif algo == "h2pc"
        R"dagout = h2pc(dataset)"  ##R obj, dag 
    else
        error("Method $algo is not implemented. Select hc, tabu, mmhc, or h2pc.")
    end
    println("done.")

    print("Scoring net structure ($score)...")
    if score == "aic"
        bns = rcopy(R"bnscore = score(dagout, data=dataset, type='aic')")    ##score dag 
    elseif score == "bic"
        bns = rcopy(R"bnscore = score(dagout, data=dataset, type='bic')")    ##score dag 
    else
        error("Select aic or bic for scoring the network.")
    end
    
    bns = round(bns, digits=2)
    
    println("done.")
    println("Network score ($score): $bns")
    
    print("Learning parameters ($method)...")
    if method == "bayes"
        fdag = rcopy(R"fitted_dag = bn.fit(dagout, dataset, method = 'bayes')") # cond.prob. tables 
    elseif method == "mle"
        fdag = rcopy(R"fitted_dag = bn.fit(dagout, dataset, method = 'mle')") # cond.prob. tables
    else
        error("Please select bayes or mle for the parameter-learning method.")
    end
    println("done.")

    asx = rcopy(R"arc.strength(dagout, data=dataset, criterion='x2')")
    asb = rcopy(R"arc.strength(dagout, data=dataset, criterion='bic')")
    asa = rename!(asx, :strength => :X2)
    asb = rename!(asb, :strength => :Network_score_change)

    print("Chi-square probability for arcs among nodes... ")
    println(asx)
    println()
    
    print("Change in current network score [$score score = $bns] if arc removed... ")
    println(asb)
    println()
    
    R"plot(dagout)"


    R"library(gRain)"
    
    R"g <- igraph.from.graphNEL(as.graphNEL(fitted_dag))" # edges for dag
    R"betweenness <- igraph::betweenness(g, directed = FALSE)" # NODE betweeness scores [int]
    R"betweenness_std <- betweenness / ((vcount(g) - 1) * (vcount(g) - 2) / 2)" #graph normalization
    R"node_betweenness <- data.frame(betweenness = betweenness, betweenness_std = betweenness_std) %>% tibble::rownames_to_column()" #df
#    R"betwe_DF = node_betweenness %>% arrange(-betweenness) %>% .[1:10, ]" #not needed?
    R"edge_betweenness <- igraph::edge_betweenness(g, directed = FALSE)"  #EDGE betweeness
    R"betw_edges = data.frame(edge = attr(E(g), 'vnames'), betweenness = edge_betweenness) %>% tibble::rownames_to_column()" #df
    R"betw_edges_df = data.frame(str_split_fixed(betw_edges$edge, '\\|', 2))" #split edge descriptor
    R"PanelAjunction_AN2 = compile(as.grain(fitted_dag), propagate = T)" #gRain make junction tree

    R"for (i in 1:nrow(betw_edges_df)){
          betw_edges_df[i,3] =  (((querygrain(setFinding(PanelAjunction_AN2, nodes = betw_edges_df[i,1], states = '0'), nodes = betw_edges_df[i,2])[[1]][[1]]/
          querygrain(setFinding(PanelAjunction_AN2, nodes = betw_edges_df[i,1], states = '1'), nodes = betw_edges_df[i,2])[[1]][[1]])+
          (querygrain(setFinding(PanelAjunction_AN2, nodes = betw_edges_df[i,2], states = '0'), nodes = betw_edges_df[i,1])[[1]][1]/
          querygrain(setFinding(PanelAjunction_AN2, nodes = betw_edges_df[i,2], states = '1'), nodes = betw_edges_df[i,1])[[1]][1]))/2)
          }"  #Iterate junction tree and get cond. probs for each edge (as V3)

    R"betw_edges_df_CondProbs = cbind(betw_edges,betw_edges_df)" #combine edges
    R"betw_edges_df_CondProbs =betw_edges_df_CondProbs[,c(1,2,3,6,4,5)]"
    R"names(betw_edges_df_CondProbs)[4:6] = c('Fold_Ratio','from','to')" #rearrange and rename
    R"betw_edges_df_CondProbs$scaled_Bet = scale(betw_edges_df_CondProbs$betweenness)" # centered,norm by std (scale func)
    R"betw_edges_df_CondProbs$Fold_Ratio = ifelse(betw_edges_df_CondProbs$Fold_Ratio=='Inf',0,betw_edges_df_CondProbs$Fold_Ratio)" #Inf=0
    R"betw_edges_df_CondProbs$scaled_Fold = scale(betw_edges_df_CondProbs$Fold_Ratio)"
    R"betw_edges_df_CondProbs$scaled_Fold = ifelse(betw_edges_df_CondProbs$scaled_Fold=='NaN',1,betw_edges_df_CondProbs$Fold_Ratio)" #NaN=1 
    R"betw_edges_df_CondProbs$Scaled_Weighted = betw_edges_df_CondProbs$scaled_Bet * betw_edges_df_CondProbs$scaled_Fold"
    R"betw_edges_df_CondProbs$Weighted_A = betw_edges_df_CondProbs$betweenness * betw_edges_df_CondProbs$Fold_Ratio"
    R"Weighted = betw_edges_df_CondProbs$Weighted"
    R"betw_edges_df_CondProbs$Weighted2 = betw_edges_df_CondProbs$Weighted_A*0.0003" #small number multiplier
    R"write.csv(betw_edges_df_CondProbs,'Edge_Betweenes_Calcs.csv',row.names = F)"

    R"data <- toVisNetworkData(g)"  #g from igraph
    R"visIgraph(g, layout = 'layout_nicely', physics = FALSE, smooth = TRUE)"
    R"df_edges = data$edges"  #weighted edges
    R"df_edges2 = merge(df_edges,betw_edges_df_CondProbs,by.x=c('from','to'),by.y = c('from','to'))"
    R"df_edges2$Weighted = as.numeric(df_edges2$Scaled_Weighted)"
    R"df_edges2$abs_weig = abs(df_edges2$Weighted)"
    R"df_edges2_final = df_edges2[,c('from','to','abs_weig')]"  #rearranged df
    R"colnames(df_edges2_final) = c('from','to','scaled_arc_weight')"
    zz = rcopy(R"df_edges2_final")
println(zz)
    R"df_edges2_final$value = ifelse(((log(((df_edges2_final$value)^2),2))+10)<0,0.1,((log(((df_edges2_final$value)^2),2))+10))" #scaling?!
    R"rbPal <- colorRampPalette(c('lightblue','darkblue'))"
    R"df_edges2_final$color <- rbPal(10)[ as.numeric(cut(df_edges2_final$value, breaks = 10))]"
    R"df_edges2_final$highlight <- 'red' "
#    dataout = rcopy(R"df_edges2")

    #R"write.csv(nod_net2,'Node_Betweeness_Calcs.csv',row.names = F)"
    R"nod_net = as.data.frame(data$nodes)"        #orig nodes
    R"betweenness2 = as.data.frame(betweenness)"  #orig betweeness
    R"nod_net <- nod_net %>% mutate(font.size = ( log( ( (betweenness) + 1.0001), 10 ) ) + 25 )" #node font size
    R"nod_net = cbind(nod_net,betweenness)"  #add betweeness
    R"nod_net <- nod_net %>% mutate(size = betweenness * 0.25)" #scale size of node prop. to NODE betweeness
    R"nod_net <- nod_net %>% mutate(color = 'orange')"  #add node color

    #Make an html object.
    R"network = visNetwork(nod_net, edges = df_edges2_final,idToLabel=T,
                     width = '100%', height = 900) %>%
              visOptions(highlightNearest = list(enabled =TRUE, degree = 0,hover=T,labelOnly=T),
                   nodesIdSelection = TRUE) %>%
              visEdges(smooth = list(enabled = TRUE, type = 'diagonalCross',roundness = 0.1),
                   physics = F)%>%
              visNodes(color = list(background = 'lightblue',
                   border = 'darkblue',
                   highlight = 'yellow'))%>%
              visInteraction(navigationButtons = TRUE,
                   dragNodes = TRUE,
                 dragView = FALSE, zoomView = FALSE)%>%
              visIgraphLayout(layout = 'layout_with_fr',type = 'full')"
    R"visSave(network, file = 'BN.html', background = '#F5F4F4')"

# Make a D3 network    
#    R"df_edges2_final_B = df_edges2_final[,c(1:3)]"
#    R"df_edges2_final_C = data.frame(
#         from = match(df_edges2_final_B$from,nod_net$id),
#         to = match(df_edges2_final_B$to,nod_net$id),
#         value = df_edges2_final_B$value)"
#    R"nod_net$label = nod_net$id"
#    R"nod_net$id = 1:nrow(nod_net)"
#    R"nod_netC <- mutate(nod_net, id = id - 1)"
#    R"nod_netC = nod_netC[,c(1,3)]"
#    R"df_edges2_final_D <- mutate(df_edges2_final_C, from = from - 1, to = to - 1)"

#    R"nodes_d3 <- mutate(nod_net, id = id - 1)"
#    R"edges_d3 <- mutate(edges, from = from - 1, to = to - 1)"

    println("Run completed --  bnlearn net written to BN.html")

    return fdag, asa, asb

end
