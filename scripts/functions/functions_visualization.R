#### Title: Functions for data visualization 
#### Author: María Bueno Álvez
#### Description: script collecting functions to plot data
#### Last edited : 12/08/2024

# Visualization packages
library(ggrepel)
library(tidytext)
library(embed)
library(ggbeeswarm)
library(patchwork)
library(ggsci)
library(ggplotify)
library(pheatmap)
library(ggridges)
library(viridis)
library(GGally)
library(ggrain)

# Function to generate PCA
do_pca <- function(data,
                   meta = NULL,
                   variable = NULL,
                   wide = T,
                   impute = T,
                   plots = F) {
  if (wide) {
    data_w <- data
  } else {
    data_w <-
      data |>
      select(Sample, Assay, NPX) |>
      pivot_wider(values_from = NPX,
                  names_from = Assay)
  }
  
  if (impute) {
    pca_rec <-
      recipe( ~ ., data = data_w) %>%
      update_role(Sample, new_role = "id")  |>
      step_normalize(all_predictors()) |>
      step_impute_knn(all_predictors()) |>
      step_pca(all_predictors())
    
    pca_prep <- prep(pca_rec)
    
    tidied_pca <- tidy(pca_prep, 3)
    
  } else {
    pca_rec <-
      recipe( ~ ., data = data_w) %>%
      update_role(Sample, new_role = "id")  |>
      step_normalize(all_predictors()) |>
      step_pca(all_predictors())
    
    pca_prep <- prep(pca_rec)
    
    tidied_pca <- tidy(pca_prep, 2)
  }
  loadings_data <-
    tidied_pca |>
    rename(Assay = terms,
           Value = value,
           PC = component)
  
  pca_res <-  juice(pca_prep)
  
  if (plots) {
    # Loadings plot
    loadings_plot <-
      tidied_pca %>%
      filter(component %in% paste0("PC", 1:4)) %>%
      group_by(component) %>%
      top_n(8, abs(value)) %>%
      ungroup() %>%
      mutate(terms = reorder_within(terms, abs(value), component)) %>%
      ggplot(aes(abs(value), terms, fill = value > 0)) +
      geom_col() +
      facet_wrap( ~ component, scales = "free_y") +
      scale_y_reordered() +
      labs(x = "Absolute value of contribution",
           y = NULL, fill = "Positive?") +
      theme_hpa()
    
    # PCA plot
    pca_plot <-
      pca_res %>%
      left_join(meta, by = "Sample") %>%
      ggplot(aes(PC1, PC2)) +
      geom_point(aes(color = !!sym(variable)), alpha = 0.7, size = 2) +
      labs(color = NULL) +
      theme_hpa() +
      labs(color = variable)
    
    return(
      list(
        "pca_res" = pca_res,
        "loadings" = loadings_data,
        "pca_plot" = pca_plot,
        "loadings_plot" = loadings_plot
      )
    )
  } else {
    return(list("pca_res" = pca_res,
                "loadings" = loadings_data))
  }
  
}

# Function to generate UMAP
do_umap <- function(data,
                    meta = NULL,
                    variable = NULL,
                    wide = T,
                    impute = T,
                    plots = F,
                    n_neighbors = 15) {
  if (wide) {
    data_w <- data
  } else {
    data_w <-
      data |>
      select(Sample, Assay, NPX) |>
      pivot_wider(values_from = NPX,
                  names_from = Assay)
  }
  
  if (impute) {
    umap_rec <-
      recipe( ~ ., data = data_w) %>%
      update_role(Sample, new_role = "id")  |>
      step_normalize(all_predictors()) |>
      step_impute_knn(all_predictors()) |>
      step_umap(all_predictors(), neighbors = n_neighbors)
    
    umap_prep <- prep(umap_rec)
    
  } else {
    umap_rec <-
      recipe( ~ ., data = data_w) %>%
      update_role(Sample, new_role = "id")  |>
      step_normalize(all_predictors()) |>
      step_umap(all_predictors(), neighbors = n_neighbors)
    
    umap_prep <- prep(umap_rec)
    
  }
  
  umap_res <-  juice(umap_prep)
  
  if (plots) {
    # Loadings plot
    umap_plot <-
      umap_res |>
      left_join(meta, by = "Sample") |>
      ggplot(aes(UMAP1, UMAP2, color = !!sym(variable))) +
      geom_point(alpha = 0.7, size = 2) +
      theme_hpa()
    
    return(list("umap_res" = umap_res,
                "umap_plot" = umap_plot))
  } else {
    return(umap_res)
  }
  
}

# Function to generate volcano plot from differential expression results                                                                                                                                  
plot_volcano <- function(de_results) {
  
  labels <- 
    de_results |> 
    top_n(n = 10, wt = -log10(adj.P.Val)) 
  
  volcano_plot <- 
    de_results |> 
    ggplot(aes(x = logFC, y = -log10(adj.P.Val), color = sig, label = Assay)) +
    geom_point(size = 1, alpha = 0.4, show.legend = F) + 
    geom_text_repel(data = labels, size = 2, show.legend = F) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = "darkgrey") +
    geom_vline(xintercept = -0.5, linetype = 'dashed', color = "darkgrey") +
    geom_vline(xintercept = 0.5, linetype = 'dashed', color = "darkgrey") +
    scale_color_manual(values = pal_de) +
    theme_hpa() +   
    theme(axis.text = element_text(size = 8),
          legend.position = "top") 
  
  return(volcano_plot)
}

# Function to generate a beesarm plot in MEDECA/ALLVOS for selected proteins
plot_beeswarm_cohorts <- function(proteins) {
  data |> 
    filter(Assay %in% proteins) |> 
    left_join(metadata |> 
                select(Sample, Cancer, Cohort), by = "Sample") |> 
    mutate(Cohort = factor(Cohort, levels = c("MEDECA", "ALLVOS")),
           Assay = factor(Assay, levels = top_proteins)) |> 
    ggplot(aes(Cancer, NPX, fill = Cancer, color = Cancer)) +
    geom_quasirandom(size = 0.5) +
    geom_boxplot(alpha = 0.5, outlier.color = NA) +
    stat_summary(fun = "median",
                 geom = "crossbar", 
                 width = 0.2,
                 colour = "grey20",
                 show.legend = F) +
    facet_grid(Assay~Cohort, scales = "free") +
    scale_color_manual(values = pal_cancer) +
    scale_fill_manual(values = pal_cancer) +
    theme_hpa()
  
}

# Function to generate a raincloud plot in MEDECA/ALLVOS for selected proteins
plot_raincloud_cohorts <- function(proteins) {
  data |> 
    filter(Assay %in% proteins) |> 
    left_join(metadata |> 
                select(Sample, Cancer, Cohort), by = "Sample") |> 
    mutate(Cohort = factor(Cohort, levels = c("MEDECA", "ALLVOS")),
           Assay = factor(Assay, levels = top_proteins)) |> 
    ggplot(aes(Cancer, NPX, fill = Cancer, color = Cancer)) +
    geom_rain(alpha = 0.5) +
    facet_grid(Assay~Cohort, scales = "free") +
    scale_color_manual(values = pal_cancer) +
    scale_fill_manual(values = pal_cancer) +
    theme_hpa()
  
}
