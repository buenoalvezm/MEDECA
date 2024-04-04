
# Utility funcitons

# Save results 
savepath <- 
  function(savename) { 
    result_folder <- paste0("results/", Sys.Date())
    dir.create(result_folder, showWarnings = FALSE)
    
    savename <-
      paste0(result_folder, "/", savename)
    
    
    return(savename)
    
  }

savepath_data <- 
  function(folder, savename) { 
    result_folder <- paste0("data/processed/", folder)
    dir.create(result_folder, showWarnings = FALSE)
    
    savename <-
      paste0(result_folder, "/", savename)
    
    
    return(savename)
    
  }


# PCA calculation
do_pca <- 
  function(wide_data, npcs = NULL) {
    suppressMessages(require(tidyverse))
    suppressMessages(require(pcaMethods))
    
    if(is.null(npcs)) {
      npcs <- min(dim(wide_data))
    }
    
    wide_data %>% 
      t() %>% 
      pca(nPcs = npcs) 
  }

get_pca_scores <- 
  function(pca_res, 
           use_R2cum_PCselection = F,
           R2cum_lim = 0.8,
           use_sDev_PCselection = F) {
    suppressMessages(require(tidyverse))
    
    pc_lim <- ncol(pca_res@scores)
    
    if(use_R2cum_PCselection) {
      pc_lim <- 
        which(pca_res@R2cum > R2cum_lim)[1]
    } else if(use_sDev_PCselection) {
      pc_lim <- 
        rev(which(pca_res@sDev >= 1))[1]
    }
    
    pca_res %>% 
      scores() %>% 
      {.[,1:pc_lim]} 
  }

# UMAP calculation
do_umap <- 
  function(wide_data, 
           seed = 42, 
           n_neighbors = 15,
           n_components = 2, 
           ...) {
    suppressMessages(require(tidyverse))
    suppressMessages(require(uwot))
    
    set.seed(seed)
    
    umap_res <- 
      umap(wide_data, 
           n_neighbors = n_neighbors,
           n_components = n_components,
           ...) %>% 
      as.data.frame()
    
    rownames(umap_res) <- rownames(wide_data)
    colnames(umap_res) <- paste0("UMAP", 1:ncol(umap_res))
    
    umap_res
  }

# Perform overrepresentation analyses
perform_ORA <-
  function(gene_associations,
           database,
           universe,
           n_clusters = 5,
           minGSSize = 10,
           maxGSSize = Inf,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.2) {
    require(clusterProfiler)
    require(multidplyr)
    
    if(n_clusters != 1) {
      worker_cluster <- new_cluster(n = n_clusters)
      cluster_library(worker_cluster, c("dplyr",
                                        "tidyverse"))
      cluster_copy(worker_cluster, c("enricher",
                                     "universe",
                                     "database",
                                     "minGSSize",
                                     "maxGSSize",
                                     "pvalueCutoff",
                                     "qvalueCutoff" ))
      
      pre_out <- 
        gene_associations %>%
        group_by(partition) %>%
        partition(worker_cluster) 
    } else {
      pre_out <- 
        gene_associations %>%
        group_by(partition)
    }
    
    outdata <-
      pre_out %>% 
      do({
        g_data <- .
        pull(g_data, gene) %>%
          enricher(universe = universe,
                   TERM2GENE = database, 
                   minGSSize = minGSSize,
                   maxGSSize = maxGSSize,
                   pvalueCutoff = pvalueCutoff,
                   qvalueCutoff = qvalueCutoff) %>%
          as_tibble()
      }) %>%
      ungroup() %>%
      collect()
    
    if(n_clusters != 1) rm(worker_cluster)
    outdata
  }

# Missing value imputation
impute_values <- 
  function(data, ID, wide_data = F) {
    
    if(wide_data == F) {
      data_wide <- 
        data %>% 
        select(ID, Assay, NPX) %>% 
        spread(Assay,NPX) 
      
    } else {
      data_wide <- 
        data
    }
    
    data_imputed <- 
      data_wide %>% 
      column_to_rownames(ID) %>% 
      as.matrix() %>% 
      t() %>% 
      impute.knn() 
    
    final_data <- 
      data_imputed$data %>% 
      t() %>% 
      as_tibble(rownames = ID)
    
    return(final_data)
    
  }
