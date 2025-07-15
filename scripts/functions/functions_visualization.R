#### Title: Functions for data visualization 
#### Author: María Bueno Álvez
#### Description: script collecting functions to plot data
#### Last edited : 15/07/2025

# Visualization packages
library(ggrepel)
library(tidytext)
library(embed)
library(ggbeeswarm)
library(patchwork)
library(ggplotify)
library(pheatmap)
library(ggridges)
library(ggforestplot)


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
    
    # Extract proportion of variance explained
    pca_variance <- 
      tidy(pca_prep, number = 3, type = "variance")|> 
      filter(terms  == "percent variance") 
    
    tidied_pca <- tidy(pca_prep, 3)
    
  } else {
    pca_rec <-
      recipe( ~ ., data = data_w) %>%
      update_role(Sample, new_role = "id")  |>
      step_normalize(all_predictors()) |>
      step_pca(all_predictors())
    
    pca_prep <- prep(pca_rec)
    
    pca_variance <- 
      tidy(pca_prep, number = 3, type = "variance")|> 
      filter(terms  == "percent variance") 
    
    tidied_pca <- tidy(pca_prep, 2)
  }
  loadings_data <-
    tidied_pca |>
    dplyr::rename(Assay = terms,
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
        "pca_variance" = pca_variance,
        "loadings_plot" = loadings_plot
      )
    )
  } else {
    return(list("pca_res" = pca_res,
                "pca_variance" = pca_variance,
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
plot_beeswarm_cohorts <- function(proteins,
                                  data,
                                  metadata,
                                  variable,
                                  palette,
                                  angled = F) {
  data |> 
    filter(Assay %in% proteins,
           Sample %in% metadata$Sample) |> 
    left_join(metadata |> 
                select(Sample, !!sym(variable), Cohort), by = "Sample") |> 
    mutate(Cohort = factor(Cohort, levels = c("MEDECA", "ALLVOS")),
           Assay = factor(Assay, levels = proteins)) |> 
    ggplot(aes(!!sym(variable), NPX, fill = !!sym(variable), color = !!sym(variable))) +
    geom_quasirandom(size = 0.5) +
    geom_boxplot(alpha = 0.5, outlier.color = NA) +
    stat_summary(fun = "median",
                 geom = "crossbar", 
                 width = 0.2,
                 colour = "grey20",
                 show.legend = F) +
    facet_grid(Assay~Cohort, scales = "free") +
    scale_color_manual(values = palette) +
    scale_fill_manual(values = palette) +
    theme_hpa(angled = angled)
  
}

plot_beeswarm_pancancer <- function(protein) {
  
  cancer <- 
    pancancer_markers |> 
    filter(Protein == protein) |> 
    pull(Cancer)
  
  data |> 
    filter(Sample %in% c(selected_samples$Sample, 
                         medeca_healthy$Sample),
           Assay == protein) |> 
    left_join(selected_samples, by = "Sample") |>
    mutate(Cancer = ifelse(is.na(Cancer), "No diagnosis", Cancer),
           Cancer = factor(Cancer, levels = cancer_groups)) |> 
    ggplot(aes(Cancer, NPX, color = Cancer, fill = Cancer)) +
    geom_quasirandom(alpha = 0.8, size = 1) +
    geom_boxplot(alpha = 0.5, outlier.color = NA) +
    facet_wrap(~Assay, scales = "free_y", nrow = 1) +
    scale_color_manual(values = pal_cancers) +
    scale_fill_manual(values = pal_cancers) +
    theme_hpa(angled = T,
              axis_x = F) +
    ggtitle(paste0(protein, " - ", cancer))
}


plot_beeswarm <- function(data,
                          metadata,
                          proteins,
                          variable,
                          palette,
                          n_rows = 1) {

  dat <- 
    data |> 
    filter(Assay %in% proteins) |> 
    mutate(Assay = factor(Assay, levels = proteins)) |> 
    left_join(metadata, by = "Sample") 
  
  beeswarm_plot <- 
    dat |> 
    ggplot(aes(!!sym(variable), NPX, color = !!sym(variable), fill = !!sym(variable))) +
    geom_quasirandom(alpha = 0.8, size = 1) +
    geom_boxplot(alpha = 0.5, outlier.color = NA) +
    stat_summary(fun = "median",
                 geom = "crossbar", 
                 width = 0.2,
                 colour = "grey20",
                 show.legend = F) +
    scale_color_manual(values = palette) +
    scale_fill_manual(values = palette) +
    theme_hpa(angled = T,
              axis_x = F) 
  
  if(length(proteins) == 1) {
    final_plot <- 
      beeswarm_plot +
      ggtitle(proteins)
  } else {
    final_plot <- 
      beeswarm_plot +
      facet_wrap(~Assay, scales = "free_y", nrow = n_rows) 
  }
  
  return(final_plot)
}

plot_heatmap <-
  function(data,
           proteins,
           metadata,
           variable) {
    ann <-
      meta_medeca |>
      select(Sample,!!sym(variable)) |>
      column_to_rownames("Sample")
    
    data |>
      filter(Sample %in% metadata$Sample,
             Assay %in% proteins) |>
      group_by(Assay) |> 
      mutate(NPX = scales::rescale(NPX, to = c(0:1))) |> 
      pivot_wider(names_from = "Assay",
                  values_from = "NPX") |>
      column_to_rownames("Sample") |>
      #scale() |>
      t() |>
      pheatmap(
        show_colnames = F,
        clustering_method = "ward.D2",
        color = pal_heatmap(100),
        annotation_col = ann
      )
  }

plot_lollipop <- function(protein_list) {
  
  medeca_de_prots <- 
    medeca_de |> 
    filter(sig != "not significant")
  
  protein_list |>
    mutate(Pan_cancer_marker = ifelse(term %in% pancancer_markers$Protein, "Yes", "No"),
           DE = ifelse(term %in% medeca_de_prots$Assay, "Yes", "No")) |>
    ggplot(aes(fct_reorder(term, abs(estimate)), abs(estimate), color = DE)) +
    geom_segment(aes(
      x = fct_reorder(term, abs(estimate)),
      xend = fct_reorder(term, abs(estimate)),
      y = 0,
      yend = abs(estimate)
    )) +
    geom_point() +
    coord_flip() +
    theme_hpa(angled = T) +
    xlab("") +
    ylab("Model estimate") +
    scale_color_manual(values = c("grey90", "grey30")) +
    theme(
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9),
      axis.ticks.y = element_blank(),
      legend.position = "top"
    )
}

plot_densities <- function(predictions,
                           test_data,
                           metadata) {
  
  densiy_data <- 
    predictions|> 
    mutate(Sample = test_data$Sample) |>
    left_join(metadata |> 
                select(Sample, Disease_type), by = "Sample") |> 
    mutate(Disease_type = factor(Disease_type, levels = c("Other","No diagnosis", "Infectious", "Autoimmune", "Inflammatory", "Cancer"))) 
  
  densiy_data |> 
    ggplot(aes(.pred_Case, fill = Disease_type,  color = Disease_type)) +
    geom_density(alpha = 0.7, show.legend = F) +
    geom_vline(data = densiy_data |> 
                 group_by(Disease_type) |> 
                 summarise(mean_pred = mean(.pred_Case, na.rm = TRUE)), 
               aes(xintercept = mean_pred, color = Disease_type), 
               linetype = "dashed", size = 1, show.legend = F) +
    scale_fill_manual(values = c(pal_controls, "Cancer" = "#E16C54")) +   
    scale_color_manual(values = c(pal_controls, "Cancer" = "#E16C54")) +   
    facet_wrap(~Disease_type, ncol = 1) +
    theme_hpa() +
    theme(legend.position = "top") +
    xlab("Probability cancer") 
  
}

plot_roc_curve <- function(roc,
                           auc) {
  roc |> 
    ggplot(aes(
      x = 1 - specificity,
      y = sensitivity
    )) +
    geom_path(size = 1, color = "grey30") +
    geom_segment(aes(
      x = 0,
      y = 0,
      xend = 1,
      yend = 1
    ),
    colour = 'grey',
    linetype = 'dotdash') +
    scale_x_continuous(breaks = c(0, 1)) +
    scale_y_continuous(breaks = c(0, 1)) +
    theme_hpa() +
    coord_fixed() +
    theme(legend.position = "top") +
    ggtitle(paste0("AUC = ", round(auc, 2)))
}

plot_roc_tiles <- function(predictions, 
                           test_data,
                           metadata) {
  controls <- 
    metadata |> 
    distinct(Disease_type) |> 
    filter(!Disease_type %in% c("Cancer", "Other")) |> 
    pull()
  
  auc_controls <- 
    map_df(controls, function(control) {
      predictions |> 
        mutate(Sample = test_data$Sample) |>
        left_join(metadata |> 
                    select(Sample, Disease_type), by = "Sample") |> 
        filter(Disease_type %in% c(control, "Cancer")) |> 
        roc_auc(truth = Class, .pred_Case, event_level = "second") |> 
        mutate(Control = control)
      
    })
  
  roc_controls <- 
    map_df(controls, function(control) {
      predictions |> 
        mutate(Sample = test_data$Sample) |>
        left_join(metadata |> 
                    select(Sample, Disease_type), by = "Sample") |> 
        filter(Disease_type %in% c(control, "Cancer")) |> 
        roc_curve(truth = Class, .pred_Case, event_level = "second") |> 
        mutate(Control = control)
      
    })
  
  tiles <- tibble(
    Control = auc_controls$Control,
    AUC = round(auc_controls$.estimate, 2)
  ) |> 
    mutate(Control = factor(Control, levels = rev(c("No diagnosis", "Infectious", "Autoimmune", "Inflammatory")))) |> 
    arrange(Control) |> 
    mutate(specificity = rep(0.2, 4),
           sensitivity = c(0.1, 0.2, 0.3, 0.4))
  # Modify the ROC plot to include facets for each disease group
  roc_controls |> 
    ggplot(aes(
      x = 1 - specificity,
      y = sensitivity,
      group = Control
    )) +
    geom_path(aes(color = Control), size = 1) +
    geom_segment(
      aes(x = 0, y = 0, xend = 1, yend = 1),
      colour = 'grey',
      linetype = 'dotdash'
    ) +
    geom_tile(
      width = 0.2, 
      height = 0.1,
      data = tiles,
      aes(fill = Control),
      show.legend = FALSE
    ) +
    geom_text(
      data = tiles,
      aes(label = AUC),
      size = 3,
      color = "white",
      show.legend = FALSE
    ) +
    geom_text(
      label = "AUC",
      x = 0.8,
      y = 0.5,
      size = 3,
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = pal_controls) + # Add colors for each disease
    scale_fill_manual(values = pal_controls) +
    scale_x_continuous(breaks = c(0, 1)) +
    scale_y_continuous(breaks = c(0, 1)) +
    theme_hpa() +
    coord_fixed() +
    theme(legend.position = "top")
}

plot_ora <- function(enrichment,
                     protein_list,
                     pval_lim = 0.05,
                     ncateg = 10,
                     fontsize = 10) {
  
  # From gene name to ENTREZID
  protein_conversion <- clusterProfiler::bitr(protein_list,
                                              fromType = "SYMBOL",
                                              toType = "ENTREZID",
                                              OrgDb = org.Hs.eg.db::org.Hs.eg.db)
  
  protein_list <- protein_conversion |> dplyr::pull(ENTREZID) |> unique()
  
  # Visualize results
  dotplot <- clusterProfiler::dotplot(enrichment,
                                      showCategory = ncateg,
                                      font.size = fontsize)
  
  barplot <- barplot(enrichment,
                     drop = TRUE,
                     showCategory = ncateg,
                     font.size = fontsize)
  
  goplot <- clusterProfiler::goplot(enrichment,
                                    showCategory = ncateg)
  
  enrichment <- clusterProfiler::setReadable(enrichment, OrgDb = org.Hs.eg.db::org.Hs.eg.db)
  cnetplot <- clusterProfiler::cnetplot(enrichment,
                                        showCategory = ncateg,
                                        categorySize = "pvalue",
                                        color.params = list(foldChange = protein_list),
                                        cex.params = list(category_label = (fontsize + 2)/12,
                                                          gene_label = (fontsize)/12))
  
  return(list("dotplot" = dotplot,
              "barplot" = barplot,
              "goplot" = goplot,
              "cnetplot" = cnetplot))
}

plot_gsea <- function(enrichment,
                      de_results,
                      pval_lim = 0.05,
                      ncateg = 10,
                      fontsize = 10) {
  
  # Prepare sorted_protein_list
  protein_list <- stats::setNames(de_results$logFC,
                                  de_results$Assay)
  sorted_protein_list <- sort(protein_list, decreasing = TRUE)
  
  # Visualize results
  dotplot <- clusterProfiler::dotplot(enrichment,
                                      showCategory = ncateg,
                                      font.size = fontsize,
                                      split=".sign") +
    ggplot2::facet_grid(.~.sign)
  
  ridgeplot <- clusterProfiler::ridgeplot(enrichment, showCategory = ncateg) +
    ggplot2::labs(x = "enrichment distribution") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = fontsize),
                   axis.text.y = ggplot2::element_text(size = fontsize),
                   text = ggplot2::element_text(size = fontsize))
  
  gseaplot <- clusterProfiler::gseaplot(enrichment,
                                        by = "all",
                                        title = enrichment$Description[1],
                                        geneSetID = 1)
  
  enrichment <- clusterProfiler::setReadable(enrichment, OrgDb = org.Hs.eg.db::org.Hs.eg.db)
  cnetplot <- clusterProfiler::cnetplot(enrichment,
                                        showCategory = ncateg,
                                        categorySize = "pvalue",
                                        color.params = list(foldChange = protein_list),
                                        cex.params = list(category_label = (fontsize + 2)/12,
                                                          gene_label = (fontsize)/12))
  
  return(list("enrichment" = enrichment,
              "dotplot" = dotplot,
              "cnetplot" = cnetplot,
              "ridgeplot" = ridgeplot,
              "gseaplot" = gseaplot))
}
