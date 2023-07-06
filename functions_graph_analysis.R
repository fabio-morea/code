# functions used in several notebooks
# 
# ############# Calculate empirical value of mixing parameter MU
# 


empirical_mu <- function(G) {
    gdf<-as_long_data_frame(G) 
    gdf$inter_comm <- (gdf$from_community != gdf$to_community)
    inter_community_links <- sum(gdf$w[ gdf$inter_comm == TRUE]) 
    mu = sum(inter_community_links) / sum(gdf$w)
    return(mu)
}


# ################# community detection, 
# retunrs a community object 
# including the algorithm 
# 
find_communities <- function(g, method, verbose = FALSE) {
    if (method == "LV"){
        comms <- cluster_louvain(g, resolution = 1.0)
    } else if (method == "FG"){
        comms <- fastgreedy.community(g)
    } else if (method == "IM"){
        comms <- infomap.community(g)
    } else if (method == "LP"){
        comms <- label.propagation.community(g)
    } else if (method == "ML"){
        comms <- multilevel.community(g)
    } else if (method == "WT"){
        comms <- walktrap.community(g)
    } else if (method == "LE"){
        comms <- leading.eigenvector.community(g)    
    } else {
        print("No valid method")
        stop
    }
    comms$algorithm = method
    
    if (verbose == TRUE) {print(paste("Community detection with ", method, "completed."))}
    return(comms)
}


find_communities_N_times <- function(g,
                                     n_trials,
                                     methods,
                                     filename_summary = '',
                                     filename_membership = '') {
    results_n_trials <- data.frame()
    membership_matrix <- c()
 
    for (i in 1:n_trials) {
        for (met in methods) {
            tmp_comms <- find_communities(g, method = met)
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
analyse_communities <- function(g, communities, verbose = FALSE) {
    method <- communities$algorithm
    c_membership <- communities$membership
    # modularity: need to use c_membership + 1 to handle community label 0
    mod <- round(modularity (g,  c_membership + 1), digits = 4)
    
    #number of communities
    nc <- length(table(c_membership))
    nc_builtin <- max(V(g)$community)
    nc_norm <- nc / nc_builtin
    
    mu_built_in <- round(empirical_mu(g), 4)
    
    # NMI against "Built In Communities"
    nmi = round(aricode::NMI(as.factor(V(g)$community), as.factor(c_membership)), 3)
    
    #empirical value of mu on communities found
    V(g)$community <- c_membership
    mu_emp <- round(empirical_mu(g), 4)
    
    
    if (verbose == TRUE) {
        print(paste("Muxing parameter Mu (empirical value): ", mu))
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


normalized_co_occurrence <- function(all_clusters) {
    n_trials <- ncol(all_clusters)
    x <-matrix(0,
               nrow = nrow(all_clusters),
               ncol = nrow(all_clusters))
    colnames(x) <- V(g)$name
    rownames(x) <- V(g)$name
    for (i in (1:n_trials)) {
        nclusters <- max(all_clusters[, i])
        for (k in 1:nclusters) {
            samecluster <- (which(all_clusters[, i] == k))
            nc <- length(samecluster)
            for (t in 1:nc) {
                for (j in 1:nc) {
                    x[samecluster[j], samecluster[t]] <-
                        x[samecluster[j], samecluster[t]] + 1
                }
            }
        }
    }
    x <- x / ncol(all_clusters)
    #diag(x) <- 0
    return (x)
}



#######################
#######################
####################### 
#######################   

consensus <- function(membership_matrix, gamma_lim = 0.6) {
    
    
    nco <- normalized_co_occurrence((membership_matrix))
    nodes_to_process <- nrow(nco)
    names<- V(g)$id
    cons_results <- data.frame(names)
    cons_results$membership = rep('-', vcount(g))
    
    tmp <- nco
    diag(tmp) <- 0
    max_by_row <- apply(tmp, 1, max)
    cons_results$gamma <- max_by_row
    
    community_label <- 1
    
    while (nodes_to_process >= 1)  {
        if (is.matrix(nco)) {
            next_block_ids <- array(which(nco[1, ] >= gamma_lim))
            next_block_names <- names[next_block_ids]
            
            for (node_name in next_block_names) {
                cons_results$membership[cons_results$name == node_name] <- community_label
            }
            nco <- nco[-next_block_ids, -next_block_ids]
            
        } else {
            cons_results$membership[cons_results$name == names[1]] <- community_label
            next_block_ids <- 1
        }
        
        community_label <- community_label + 1
        nodes_to_process <- nodes_to_process - length(next_block_ids)
        names<-names[-next_block_ids]
        
        
    }
    #community object
    cons_communities <- make_clusters(g, array(as.numeric(cons_results$membership)))
    cons_communities$algorithm <- paste0(method,"_cons")
    cons_communities$gamma<-cons_results$gamma
    return(cons_communities)
}




########################
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
    selected_node_id <- which(V(g)$name == selected_node_name)
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