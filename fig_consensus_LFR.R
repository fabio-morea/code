## ----setup, include=FALSE--------------------------------------------------------------------
knitr::opts_chunk$set(include = TRUE)


## ----load libraries, include=FALSE-----------------------------------------------------------
library(igraph)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(aricode) # NMI 

#install.packages("remotes")
#remotes::install_github("schochastics/netUtils")
library(netUtils)
library(VGAM)

library(devtools)  
#devtools::install_github("fabio-morea/CCD", force = TRUE)
library(CCD)
print("packages loaded")



## --------------------------------------------------------------------------------------------
make_benchmark_network <- function(n,
                                   average_degree,
                                   max_degree,
                                   min_community,
                                   max_community,
                                   mu,
                                   tau1,
                                   tau2, 
                                   max_trials=100) {
    err_check <- TRUE
    gLFR <- NULL
    trial <- 0
    while (err_check & trial < max_trials) {
        trial <- trial + 1
        tryCatch(
            gLFR <- sample_lfr(
                n = n,
                average_degree = average_degree,
                max_degree = max_degree,
                min_community = min_community,
                max_community = max_community,
                mu = mu ,
                tau1 = tau1,
                tau2 = tau2
            ),error = function(e)NULL)
        if (is.null(gLFR)) {
            print(paste("unable to generate netowork - trial", trial))
        } else {
            err_check <- FALSE
        }
    }
    if (trial > max_trials) {
        print(paste("unsuccessfully tried",max_trials,"times. Please provide a different set of params"))
        return(-1)
    } else {
        return (gLFR)
    }
    
}

testG <- make_benchmark_network(
            n = 1000,
            average_degree = 5,
            max_degree = 20,
            min_community = 5,
            max_community = 50,
            mu = 0.2,
            tau1 = 3.0,
            tau2 = 2.0
            
        )
print(testG)


## --------------------------------------------------------------------------------------------

load_benchmark_network <- function(mui, path = path, verbose = FALSE) {
  filename = paste0(path, mui, ".gml")
  g <- read_graph(filename, format = "gml")
  
  g<-as.undirected( g , mode = 'each')
  g<-igraph::simplify(g, remove.multiple = TRUE, edge.attr.comb = igraph_opt("sum"))
  
  V(g)$comm_built_in <-  V(g)$community
  V(g)$core <- coreness(g)
  V(g)$str <- strength(g)
  V(g)$name <- paste0("V" , V(g)$name)
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


## --------------------------------------------------------------------------------------------

recursive_consensus <- function(M, threshold , met, t) {
  j = 1
  
  while (j <= 10) {
    co_occ_matrix <-  CCD::normalized_co_occurrence((M))
    
    # all entries of co occurrenc ematrix below a threshold are set to zero
    co_occ_matrix[co_occ_matrix < threshold] <- 0
    
    g2 <- graph_from_adjacency_matrix(
        co_occ_matrix,
        diag = FALSE,
        weighted = TRUE,
        mode = "upper"
      )
    
    cons_communities_2 <- CCD::find_communities(g2, method = met, verbose = FALSE)
    
    
    if (length(table(E(g2)$weight)) == 1) {
      print("reached consensus")
      return(cons_communities_2)
    } else {
      print(paste("iteration", j))
    } 
    
    M  <- CCD::find_communities_repeated (g2, n_trials = t, method = met)
    
    j <- j + 1
    
    
  }
  print("reached max iterations")
  return(cons_communities_2)
}



## --------------------------------------------------------------------------------------------

path <- "./benchmark/LFR/LFR_benchmark_"
filename <- './results/results_cons_comm_detection_LFR_base.csv'
mui = seq(5, 55, 5)
n_repetitions<-1
n_trials = 50
threshold = 0.60
 
methods = c('LV', 'LD', 'LP', 'IM', 'WT')
summary_n_trials  <- data.frame()
nc_builtin <- c()



uncertain_nodes <- 0.0
uncertainty_q10 <- 0.0
uncertainty_q50 <- 0.0
uncertainty_q90 <- 0.0

for (j in 1:n_repetitions) {
    print(paste("REPETITION ", j))
    for (mu in mui) {
      print(paste("reading g ----  MU , ITERATION = ", mu, j))
      gfilename <- paste0("./benchmark/LFR/LFR_benchmark_",mu,".gml")
      g <- read_graph(gfilename, format = "gml")
      
      
      # check the main component
      components <- clusters(g)
      print("components: ")
      print(table(components$csize))
      
    V(g)$name <- paste0("V", 1:vcount(g))
    E(g)$weight <- 1.0
    E(g)$w <- 1.0
    V(g)$comm_built_in <- V(g)$membership
    V(g)$community <- V(g)$membership
    mu_built_in = round(empirical_mu(g), 3)
    nc_builtin <- max(V(g)$membership)
    print(nc_builtin)
    
    for (method_base in methods) {
      print(paste('mu: ', mu, '**** method : ', method_base))
      method <- paste0(method_base, '_ST')
      
      M <- data.frame(name = V(g)$name)
      for (i in 1:n_trials) {
        gs <-igraph::permute(g, sample(1:vcount(g), size = vcount(g), replace = FALSE))
        comms <- CCD::find_communities(gs, method = method, r=1.0, s = 5)
        # order() function to get the communities found on gs in the order they appear in g
        V(g)$community <-
          comms$membership[order(match(comms$name, V(g)$name))]
        
        nmi <-
          aricode::NMI(as.factor(V(g)$community),
                       as.factor(V(g)$comm_built_in))
        mu_emp <- round(empirical_mu(g), 3)
        nc <- length(unique(V(g)$community))
        nc_norm <- round( nc/ nc_builtin, 4)
        
        met <- method
      results_iteration <-data.frame(mu, mu_built_in, mu_emp, met, j, nmi, nc_norm, uncertain_nodes, uncertainty_q10, uncertainty_q50, uncertainty_q90)
        summary_n_trials <-
          rbind(summary_n_trials , results_iteration)
        
        comm_labeled <-
          data.frame(name = V(gs)$name,
                     memb = comms$membership)
        M <- inner_join(M, comm_labeled, by = "name")
        colnames(M) <- c("name", seq(1:i))
        
      }
      
      method <- paste0(method_base, '_consRec')
      commsRC <- recursive_consensus (M,
                                      threshold = threshold ,
                                      met = method_base,
                                      t = n_trials)
      
      V(g)$community <-
        commsRC$membership[order(match(commsRC$names, V(g)$name))]
      nmi <-
        aricode::NMI(as.factor(V(g)$community), as.factor(V(g)$comm_built_in))
      mu_emp <- round(empirical_mu(g), 3)
      nc <- length(unique(V(g)$community))
      nc_norm <- round( nc/ nc_builtin, 4)
      print(paste("comms: ",nc, nc_builtin, nc_norm))
      
      met <- method
      results_iteration <-data.frame(mu, mu_built_in, mu_emp, met, j, nmi, nc_norm, uncertain_nodes, uncertainty_q10, uncertainty_q50, uncertainty_q90)
      summary_n_trials <-
        rbind(summary_n_trials , results_iteration)
      
      
      
      method <- paste0(method_base, '_CCD_p06')
      commsCCD <- CCD::consensus_community_detection(g, p=0.6, r=1.0, group_outliers = FALSE)
      V(g)$community <- commsCCD$membership[order(match(commsCCD$name, V(g)$name))]
      uncertain_nodes <- sum(commsCCD$gamma>0)/vcount(g)
      uncertainty_q10 <- round(quantile(commsCCD$gamma[ commsCCD$gamma >0], 0.10),6)
      uncertainty_q50 <- round(quantile(commsCCD$gamma[ commsCCD$gamma >0], 0.50),6)
      uncertainty_q90 <- round(quantile(commsCCD$gamma[ commsCCD$gamma >0], 0.90),6)
      nmi <- aricode::NMI(as.factor(V(g)$community), as.factor(V(g)$comm_built_in))
      mu_emp <- round(empirical_mu(g), 3)
      nc <- length(unique(V(g)$community))
      nc_norm <- round( nc/ nc_builtin, 4)
      print(paste("comms: ",nc, nc_builtin, nc_norm))
      
      met <- method
      results_iteration <-data.frame(mu, mu_built_in, mu_emp, met, j, nmi, nc_norm, uncertain_nodes, uncertainty_q10, uncertainty_q50, uncertainty_q90)
      summary_n_trials <-rbind(summary_n_trials , results_iteration)
      print(paste(method, sum(table(commsCCD$membership) == 1), nc_norm))
      
      method <- paste0(method_base, '_CCD_p08')
      commsCCD <- CCD::consensus_community_detection(g, p=0.8, r=1.0,  group_outliers = FALSE)
      V(g)$community <- commsCCD$membership[order(match(commsCCD$name, V(g)$name))]
      
      uncertain_nodes <- sum(commsCCD$gamma>0)/vcount(g)
      uncertainty_q10 <- round(quantile(commsCCD$gamma[ commsCCD$gamma >0], 0.10),6)
      uncertainty_q50 <- round(quantile(commsCCD$gamma[ commsCCD$gamma >0], 0.50),6)
      uncertainty_q90 <- round(quantile(commsCCD$gamma[ commsCCD$gamma >0], 0.90),6)
      
      nmi <- aricode::NMI(as.factor(V(g)$community), as.factor(V(g)$comm_built_in))
      mu_emp <- round(empirical_mu(g), 3)
      nc <- length(unique(V(g)$community))
      nc_norm <- round( nc/ nc_builtin, 4)
      print(paste("comms: ",nc, nc_builtin, nc_norm))
      
      met <- method
      results_iteration <-data.frame(mu, mu_built_in, mu_emp, met, j, nmi, nc_norm, uncertain_nodes, uncertainty_q10, uncertainty_q50, uncertainty_q90)
      summary_n_trials <-rbind(summary_n_trials , results_iteration)
      print(paste(method, sum(table(commsCCD$membership) == 1), nc_norm))
      
      method <- paste0(method_base, '_CCD_p08g')
      #commsCCD <- CCD::consensus_community_detection(g, p=0.8, r=1.0, s = 5, group_outliers = TRUE
      commsCCD$membership[commsCCD$comm_size == 1] <- 0 #group single node communities
      V(g)$community <- commsCCD$membership[order(match(commsCCD$name, V(g)$name))]
      uncertain_nodes <- sum(commsCCD$gamma>0)/vcount(g)
      uncertainty_q10 <- round(quantile(commsCCD$gamma[ commsCCD$gamma >0], 0.10),6)
      uncertainty_q50 <- round(quantile(commsCCD$gamma[ commsCCD$gamma >0], 0.50),6)
      uncertainty_q90 <- round(quantile(commsCCD$gamma[ commsCCD$gamma >0], 0.90),6)
      
      nmi <- aricode::NMI(as.factor(V(g)$community), as.factor(V(g)$comm_built_in))
      mu_emp <- round(empirical_mu(g), 3)
      nc <- length(unique(V(g)$community))
      nc_norm <- round( nc/ nc_builtin, 4)
      print(paste("comms: ",nc, nc_builtin, nc_norm))
      
      met <- method
      results_iteration <-data.frame(mu, mu_built_in, mu_emp, met, j, nmi, nc_norm, uncertain_nodes, uncertainty_q10, uncertainty_q50, uncertainty_q90)
      summary_n_trials <-rbind(summary_n_trials , results_iteration)
      print(paste(method, sum(table(commsCCD$membership) == 1), nc_norm))
      print(paste("mean gamma", mean(commsCCD$gamma)))
      print(paste("NMI = ", nmi))
      uncertain_nodes <- 0.0
      uncertainty_q10 <- 0.0
      uncertainty_q50 <- 0.0
      uncertainty_q10 <- 0.0
    }
  }
    summary_n_trials  %>% write_csv(filename)
    
}

summary_n_trials  %>% write_csv(filename)

 

