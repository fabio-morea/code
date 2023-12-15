# functions used in several notebooks
# 
# 
# 

load_benchmark_network <- function(mui, path = path, verbose = FALSE) {
    filename = paste0(path, mui, ".gml")
    g <- read_graph(filename, format = "gml")
    
    g<-as.undirected( g , mode = 'each')
    g<-igraph::simplify(g, remove.multiple = TRUE, edge.attr.comb = igraph_opt("sum"))
    
    V(g)$comm_built_in <-  V(g)$community
    V(g)$core <- coreness(g)
    V(g)$str <- strength(g)
    V(g)$name <- paste0("V" , V(g)$label)
    E(g)$w <- 1.0
    E(g)$ww <- 1.0
    E(g)$weight <- 1.0
    
    # print
    if (verbose == TRUE) {
        print(paste0("Loaded benchmark network ", path, mui, ".gml"))
        print(paste("Giant component size :", length(V(g)$community)))
        print(paste("Built-in communities: ", max(V(g)$community)))
        mod <- round( modularity(g, array(V(g)$community+1) ), digits = 4)
        print(paste("Modularity of built-in communities: ", mod))
    } 
    return(g)
}

# ############# Calculate empirical value of mixing parameter MU
 
empirical_mu <- function(g) {
    gdf<-as_long_data_frame(g) 
    gdf$inter_comm <- (gdf$from_community != gdf$to_community)
    inter_community_links <- sum(gdf$weight[ gdf$inter_comm == TRUE]) 
    mu = sum(inter_community_links) / sum(gdf$weight)
    return(mu)
}


# ################# community detection, 
# retunrs a community object 
# including the algorithm 

find_communities <- function(g,method, r = c(1.0),  s = c(4), verbose = FALSE) {
    # applies selected method
    # applies resolution for LV and LD
    # undirected(g) for LV and LD
    # 
   
    method = substr(method, 1, 2)
    if (method == "LV") {
        comms <- cluster_louvain(as.undirected(g, mode = 'each'), resolution = sample(r, 1))
    } else if (method == "ML") {
        comms <- multilevel.community(as.undirected(g, mode = 'each'), resolution = sample(r, 1))
    } else if (method == "LD") {
        comms <-
            cluster_leiden(as.undirected(g, mode = 'each'), resolution_parameter = quantile(strength(g))[2] / (gorder(g) - 1))
    } else if (method == "FG") {
        comms <- fastgreedy.community(g)
    } else if (method == "IM") {
        comms <- infomap.community(g)
    } else if (method == "LP") {
        comms <- label.propagation.community(g)
    
    } else if (method == "WT") {
        comms <- walktrap.community(g, steps = sample(s, 1))
    } else if (method == "LE") {
        comms <- leading.eigenvector.community(g)
    
    } else {
        print("No valid method")
        stop
    }
    comms$algorithm = method
    
    if (verbose == TRUE) {
        print(paste("Community detection with ", method, "completed."))
    }
    return(comms)
}



find_communities_repeated <- function(g,
                                      n_trials,
                                      method = method,
                                      shuffle = TRUE,
                                      resolution = c(1.0),
                                      steps = c(4),
                                      verbose = FALSE) {
    # g must have a non null vertex attribute 'name'
    membership_table <- data.frame(name = V(g)$name)
    
    for (i in 1:n_trials) {
        if (verbose){print(i)}
        if (shuffle == TRUE) {
            gs <- igraph::permute(g, sample(1:vcount(g),size = vcount(g),replace = FALSE))
        } else {
            gs <- g
        }
        comms <-find_communities(gs, method = method, r = resolution, s = steps)
        comm_labeled <-data.frame(name = V(gs)$name, memb = comms$membership)
        membership_table <-inner_join(membership_table ,  comm_labeled, by = 'name')
        colnames(membership_table) <- c('name', seq(1:i))
    }
    membership_table <- membership_table %>% arrange(name)
    return(membership_table)
}


find_communities_N_times <- function(g,
                                     n_trials,
                                     methods,
                                     shuffle = FALSE,
                                     resolution = c(1.0),
                                     steps = c(4),
                                     filename_summary = '',
                                     filename_membership = '') {
    results_n_trials <- data.frame()
    membership_matrix <- c()
    
    
    for (i in 1:n_trials) {
        for (met in methods) {
            tmp_comms <- find_communities(g, resolution = resolution, steps = steps, method = met, shuffle = shuffle)
            V(g)$community <- tmp_comms$membership
            results_single <- analyse_communities(g, tmp_comms)
            results_n_trials = rbind(results_n_trials, results_single)
            membership_matrix <- rbind(membership_matrix, t(tmp_comms$membership))
            
            
        }
    }
    
    membership_n_trials <- data.frame(membership_matrix)
    if (filename_membership != '') {
        membership_n_trials %>% write_csv(filename_membership)
        results_n_trials %>% write_csv(filename_summary)
    }
    return(membership_n_trials)
    
    
}


# ################# analysis of communities
# returns modularity, number of communities and NMI against true labels
# communities <- cons_communities_2
analyse_communities <- function(g, communities, verbose = FALSE) {
    method <- communities$algorithm
    c_membership <- communities$membership
    # modularity: need to use c_membership + 1 to handle community label 0
    mod <- round(modularity (g,  c_membership + 1), digits = 4)
    
    #number of communities
    nc <- length(table(c_membership))
    nc_builtin <- length(table(V(g)$comm_built_in))
    nc_norm <- nc / nc_builtin
    
    mu_built_in <- round(empirical_mu(g), 4)
    
    # NMI against "Built In Communities"
    nmi = round(aricode::NMI(as.factor(V(g)$comm_built_in), as.factor(c_membership)), 3)
    
    #empirical value of mu on communities found
    V(g)$community <- c_membership
    mu_emp <- round(empirical_mu(g), 4)
     
    if (verbose == TRUE) {
        print(paste("Muxing parameter Mu (empirical value): ", mu_emp))
        print(paste("Communities found: ", nc))
        print(paste("Modularity: ", mod))
        print(paste("Normalized Mutual Information between C.D. method and built-in communities:",nmi)
        )
    }
    
    return(data.frame(method , mu_built_in, mu_emp, mod, nc, nc_norm, nmi))
}

######################
######################
######################


normalized_co_occurrence <- function(M) {
    
    names <- M$name
    M<-as.matrix(M %>% select(-name))
    n_trials <- ncol(M)
    n_nodes <- nrow(M)
    CO <-matrix(0, nrow = n_nodes,ncol =n_nodes)
    colnames(CO) <- names
    rownames(CO) <- names
    for (t in (1:n_trials)) {
        nclusters <- max(M[, t])
        for (k in 1:nclusters) {
            samecluster <- (which(M[, t] == k))
            nc <- length(samecluster)
            for (i in 1:nc) {
                for (j in 1:nc) {
                    CO[samecluster[j], samecluster[i]] <-
                        CO[samecluster[j], samecluster[i]] + 1
                }
            }
        }
    }
    # X is a matrix of x_ij counts of co-occurrence
    # X_normalized is a matrix of gamma_ij normalized coefficients
    # expressing the probability that i and j are in the same community
    # X_normalized_ij == 1.0 means that i and j have ALWAYS been classified in the same community
    # X_normalized_ij == 0.0 means that i and j have NEVER been classified in the same community
    X_normalized <- CO / n_trials
    return (X_normalized)
}



#######################
#######################
####################### 
#######################   



consensus_communities <- function(nco, gamma_lim){
    
    results <- data.frame(name = as.numeric(colnames(nco)))
    results$done <- FALSE
    results$cons_comm_label <- 0
    
    # GAMMA 
    coeffs <- nco 
    diag(coeffs)<-0.0
    results$gamma <-  round( 1 -  apply(coeffs, 1, max), 4)
    
    community_label <- 0
    nodes_to_process = nrow(results)
    while ( nodes_to_process >= 1 )  {
        community_label <- community_label + 1
        row_test <- nco[ which.max(results$done == FALSE), ]
        nodes_above_threshold <-  (row_test > gamma_lim )
        results$cons_comm_label[ nodes_above_threshold ] <- community_label
        results$done[ nodes_above_threshold] <-  TRUE
        nodes_to_process <- sum( results$done == FALSE)  
    }
    
    return(results %>% arrange(name))
}




########################


consensus_community_detection <- function(g, t, method='LV', gamma_lim, resolution=c(1.0), shuffle=TRUE) {
 
    M <- find_communities_repeated(g,
                                   n_trials=t, 
                                   method = method,
                                   shuffle = shuffle,
                                   resolution = resolution, 
                                   verbose = FALSE)
    
    nco <- normalized_co_occurrence(M)
    
    CC <- consensus_communities(nco,gamma_lim=gamma_lim)
    
    cons_communities <- make_clusters(g, array(as.numeric(CC$cons_comm_label)))
    cons_communities$gamma<-CC$gamma
    cons_communities$name<-CC$name
    return(cons_communities)
}

########################
########################

make_community_network <- function (g) {

edges_list <- g %>%
    as_long_data_frame() %>%
    select(from_community, to_community, w) %>%
    group_by(from_community, to_community) %>%    
    summarize(weight = sum(w))

gc <- graph_from_data_frame(edges_list, directed = FALSE)
V(gc)$id <- (1:vcount(gc))

comms <- data.frame(label = V(gc)$name)
comms$size <- 0
for (i in 1:length(comms$label)) {
    comms$size[i] <-
        length(V(g)$community[V(g)$community == comms$label[i]])
}
V(gc)$size <- comms$size

return(gc)
}










find_neighbours <- function(g, node_name, order) {
    selected_node_id <- which(V(g)$name == node_name)
    list_neis <- ego(g, order = order, selected_node_id)
    nei_ids <- c()
    for (nn in list_neis) {
        for (x in nn) {
            nei_ids <- append(nei_ids, x)
        }
    }
    nei_ids <- c(nei_ids, selected_node_id)
    nei_ids <- unique(nei_ids)
    return(nei_ids)
}


make_ring_of_cliques <- function(num_cliques,
                                 clique_size,
                                 add_center = TRUE,
                                 add_bridges = TRUE ) {
    f <- as.undirected(graph.empty())
    
    for (i in 1:num_cliques) {
        next_clique <- make_clique(clique_size, comm_label = paste0("C", i))
        f <- f + next_clique
    }
    
    b <- vcount(f)
    
    
    if (add_bridges) {
        f <- add_vertices(f, num_cliques)
    }
    
    
    for (j in 1:(num_cliques)) {
        b <- b + 1
        b_start <- (j-1) * clique_size +1
        b_end <- b_start + clique_size +1
        if (b_end > (clique_size * num_cliques)) {b_end <- 2}
        if (add_bridges) {
            f <- add_edges(f, c(b_start, b))
            f <- add_edges(f, c(b, b_end))
            V(f)$community[b] <- paste0("B", j)
        } else {
            f <- add_edges(f, c(b_start, b_end))
        }
        
    }
    
    if (add_center) {
        f <- add_vertices(f, 1)
        id_center <- vcount(f)
        V(f)$community[id_center] <- "A"
        for (j in 1:(num_cliques)) {
            c_start <- (j-1) * clique_size +3
            f <- add_edges(f, c(c_start , id_center))
        }
        
    }
    
    E(f)$weight <- 1.0
    
    V(g)$name <- 1:vcount(g)
    V(g)$id <-   1:vcount(g)
    
    
    return(g)
    
}