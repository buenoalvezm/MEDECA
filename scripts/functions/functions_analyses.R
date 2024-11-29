#### Title: Functions for data analyses 
#### Author: María Bueno Álvez
#### Description: script collecting functions to perform data analyses (DE & ML)
#### Last edited : 12/08/2024

# Data analyses packages
library(limma)
library(tidymodels)
library(themis)
library(eulerr)

# Function to run differential expression using limma
do_limma <-
  function(data_wide, 
           metadata,
           variable,
           case,
           correct = T) {
    
    # Filter data (cohort, variable is 0 or 1)
    dat <-
      data_wide %>% 
      inner_join(metadata %>% 
                   select(Sample, Sex, Age, !!sym(variable)), by = "Sample") %>% 
      rename(Group = !!sym(variable)) %>% 
      mutate(Group = ifelse(Group == case, "1_Case", "0_Control")) 
    
    # Design a model - add Group, and Sex, Age
    if(correct == T) {
      design <- model.matrix(~0 + as.factor(dat$Group) + as.factor(dat$Sex) + dat$Age) 
      colnames(design) <- c("control", "case",  "Sex", "Age") 
    } else {
      design <- model.matrix(~0 + as.factor(dat$Group))
      colnames(design) <- c("control", "case")
    }
    
    # Make contrast
    contrast <- makeContrasts(Diff = case - control, levels = design)
    
    dat_fit <- 
      dat %>% 
      select(-Sex, -Age, -Group)  %>% 
      column_to_rownames("Sample") %>% 
      t()
    
    fit <- lmFit(dat_fit, design = design,  method = "robust", maxit = 10000)
    
    # Apply contrast
    contrast_fit <- contrasts.fit(fit, contrast)
    
    # Apply empirical Bayes smoothing to the SE
    ebays_fit <- eBayes(contrast_fit)
    
    # Extract DE results
    DE_results <-
      topTable(ebays_fit,
               n = nrow(ebays_fit$p.value), 
               adjust.method = "fdr", 
               confint = TRUE)
    
    DE_res <- 
      DE_results %>% 
      as_tibble(rownames = "Assay") %>% 
      mutate(sig = case_when(adj.P.Val < 0.05 & logFC < -0.5 ~ "significant down",
                             adj.P.Val < 0.05 & logFC > 0.5 ~ "significant up", 
                             T ~ "not significant")) 
    
    return(DE_res)
  }


# Function to select proteins based on AUC
step_filter_auc <-  function(data, 
                             cutoff = 0.8) {

  auc_res <- 
    map_df(unique(data$Assay), function(protein) {
    protein_exp <- 
      data |> 
      filter(Assay == protein) |> 
      pull(NPX)
    
    cancer_class <- 
      data |> 
      filter(Assay == protein) |> 
      pull(Cancer)
    
    res <- pROC::roc(cancer_class, protein_exp, verbose = F)
    
    tibble(Assay = protein,
           AUC = round(as.numeric(res$auc), 2))
    
  })
  
  plot_auc <- 
    auc_res |> 
    ggplot(aes(AUC)) +
    geom_histogram() +
    geom_vline(xintercept = cutoff, color = "red") +
    theme_hpa() +
    ggtitle("AUC distribution")
  
  auc_res_cutoff <- 
    auc_res |> 
    filter(AUC > cutoff)
  
  return(list(full_list = auc_res,
              auc_distribution = plot_auc,
              final_list = auc_res_cutoff))
  
}


step_filter_t_test <- 
  function(data,
           output = "NPX",
           cutoff = 1) {
    
    ttest_results <- 
      map_df(unique(data$Assay), function(protein) {
        
        group1<-
          data %>%
          filter(Assay == protein,
                 Cancer == "Yes") %>%
          pull(NPX)
        
        group2<-
          data %>%
          filter(Assay == protein,
                 Cancer == "No") %>%
          pull(NPX)
        
        test_res <- t.test(group1 ,group2)
        p.val<- test_res$p.value
        difference <- mean(group1)-mean(group2)
        
        tibble(Assay = protein, 
               pval = p.val, 
               NPX_difference = difference)
        
      }) |> 
      mutate(pval_adj = p.adjust(pval, method = "BH")) 
    
    if(output == "NPX") {
      
      npx_res <- 
        ttest_results |> 
        select(Assay, NPX_difference)
      
      plot_npx <- 
        npx_res |> 
        ggplot(aes(NPX_difference)) +
        geom_histogram() +
        geom_vline(xintercept = cutoff, color = "red") +
        theme_hpa() +
        ggtitle("NPX difference distribution")
      
      npx_res_cutoff <- 
        npx_res |> 
        filter(NPX_difference > cutoff)
      
      return(list(full_list = npx_res,
                  npx_distribution = plot_npx,
                  final_list = npx_res_cutoff))
      
    } else if (output == "pval") {
      
      pval_res <- 
        ttest_results |> 
        select(Assay, pval_adj)
      
      plot_pval <- 
        pval_res |> 
        ggplot(aes(pval_adj)) +
        geom_histogram() +
        geom_vline(xintercept = cutoff, color = "red") +
        ggtitle("p-value distribution") +
        theme_hpa()
      
      pval_res_cutoff <- 
        pval_res |> 
        filter(pval_adj < cutoff)
      
      return(list(full_list = pval_res,
                  pval_distribution = plot_pval,
                  final_list = pval_res_cutoff))
    
    } else {
      warning("Output not recognized")
    }
    
  }


# Function to perform lasso analyses
do_glmnet <-  
  function(variable,
           case, 
           split_train, 
           split_test) {
    
    # Prepare data - make custom split for current variable
    training_dat <- 
      split_train |> 
      mutate(Class = case_when(!!sym(variable) == case ~ "Case",
                               !!sym(variable) != case ~ "Control")) |> 
      mutate(Class = factor(Class, levels = c("Control", "Case"))) |>  
      select(-!!sym(variable)) |> 
      mutate(Sample = as.character(Sample))
    
    testing_dat <- 
      split_test  |> 
      mutate(Class = case_when(!!sym(variable) == case ~ "Case",
                               !!sym(variable) != case ~ "Control")) |> 
      mutate(Class = factor(Class, levels = c("Control", "Case"))) |>  # Explicitly set levels here
      select(-!!sym(variable))|> 
      mutate(Sample = as.character(Sample))
    
    
    ml_split_custom <- make_splits(training_dat, testing_dat)
  
    
    # Define recipe
    ml_recipe <- 
      recipe(Class ~ ., data = training_dat) |> 
      #step_relevel(Class, ref_level = "Control") |> 
      update_role(Sample, new_role = "id") |> 
      step_normalize(all_numeric_predictors()) |> 
      step_nzv(all_numeric_predictors()) |> 
      step_impute_knn(all_numeric_predictors()) 
    
    # LASSO model specifications
    glmnet_specs <- 
      logistic_reg() |> 
      set_mode("classification") |> 
      set_engine("glmnet") |> 
      set_args(penalty = tune(), 
               mixture = 0) 
    
    # ML workflow
    glmnet_wflow <-
      workflow() |> 
      add_recipe(ml_recipe) |> 
      add_model(glmnet_specs) 
    
    # Define glmnet grid
    set.seed(213)
    glmnet_grid <-
      glmnet_wflow |>
      extract_parameter_set_dials() |>
      grid_latin_hypercube(size = 100)
    
    # Define the resamples (CV)
    set.seed(213)
    ml_rs <- vfold_cv(training_dat, v = 10, strata = Class)
    
    # Define the evaluation metrics (add brier)
    eval_metrics <- metric_set(roc_auc)
    
    # Define control_grid
    set.seed(213)
    ctrl <- control_grid(save_pred = TRUE, 
                         event_level = "second",
                         parallel_over = "everything") # look at extract = identity
    
    # Glmnet grid search
    set.seed(213)
    glmnet_res <-
      glmnet_wflow |>
      tune_grid(
        resamples = ml_rs,
        grid = glmnet_grid,
        control = ctrl,
        metrics = eval_metrics
      )

    plot_train <- autoplot(glmnet_res)
    
    predictions_train <- 
      glmnet_res |> 
      collect_predictions()
    
    metrics_train <- 
      glmnet_res |> 
      collect_metrics()
    
    # Select best hyperparameter
    best_glmnet <- 
      select_best(glmnet_res, metric = "roc_auc") |> 
      select(-.config)
    
    #Finalize the workflow and fit the final model
    glmnet_wflow <- 
      glmnet_wflow |>  
      finalize_workflow(best_glmnet)
    
    final_glmnet_fit <- last_fit(glmnet_wflow,
                                 ml_split_custom, 
                                 metrics = eval_metrics) 
    
    # Extract model performance
    performance <- 
      final_glmnet_fit |> 
      collect_metrics() |> 
      select(-.config, -.estimator)
    
    glmnet_auc <- 
      final_glmnet_fit |> 
      collect_metrics() |> 
      filter(.metric == "roc_auc") |> 
      pull(.estimate) |> 
      round(2)
    
    # Extract protein importance
    important_proteins <- 
      final_glmnet_fit |> 
      extract_fit_parsnip()  |> 
      vip::vi(lambda = best_glmnet$penalty, event_level = "second")  |> 
      mutate(
        Importance = abs(Importance),
        Variable = fct_reorder(Variable, Importance)
      )
    
    # Extract model predictions
    predictions <- 
      final_glmnet_fit |> 
      collect_predictions(summarize = F) 
    
    # Confusion matrix
    cm <-
      predictions |>
      mutate(pred = ifelse(.pred_Control > 0.5, "Control", "Case"),
             Class = ifelse(Class == "Control", "Control", "Case")) |>
      mutate(Class = factor(Class, levels = c("Case", "Control")),
             pred = factor(pred, levels = c("Case", "Control"))) |> 
      conf_mat(Class, pred)
    
    # ROC curve
    roc <- 
      predictions |>
      roc_curve(truth = Class, .pred_Case, event_level = "second") 
    
    
    return(list("penalty" = best_glmnet,
                "glmnet_model" = glmnet_res,
                "predictions_train" = predictions_train, 
                "performance_train" = metrics_train,
                "current_training" = training_dat, 
                "final_workflow" = glmnet_wflow,
                "final_fit" = final_glmnet_fit,
                "predictions" = predictions,
                "performance" = performance,
                "confusion_matrix" = cm,
                "roc_curve" = roc, 
                "important_proteins" = important_proteins))
  }

# Function to run protein selection and lasso models
run_ml_seed <- 
  function(seed = 213) {
    
    # Prepare split
    set.seed(seed)
    medeca_split <-
      ml_data |> 
      initial_split(prop = 0.7, strata = Cancer)
    
    medeca_train <- training(medeca_split)
    medeca_test <- testing(medeca_split)
    
    # Select proteins
    data_select_proteins <- 
      data |> 
      filter(Sample %in% medeca_train$Sample) |> 
      left_join(metadata |> 
                  select(Sample, Cancer), by = "Sample")
    
    selected_proteins_auc <- 
      step_filter_auc(data = data_select_proteins,
                      cutoff = 0.7)
    
    selected_proteins_npx <- 
      step_filter_t_test(data = data_select_proteins,
                         output = "NPX",
                         cutoff = 1)
    
    selected_protein_pval <- 
      step_filter_t_test(data = data_select_proteins,
                         output = "pval",
                         cutoff = 1E-6)
    
    y <- list("AUC > 0.7" = selected_proteins_auc$final_list$Assay, 
              "NPX diff > 1" = selected_proteins_npx$final_list$Assay,
              "adj pval < 1E-6" = selected_protein_pval$final_list$Assay)
    
    euler_plot <- plot(euler(y, shape = "circle"), quantities = TRUE, fills = list(fill = c("#449395", "#EA7C68", "#FDA638")),  alpha = 0.4) |> as.ggplot() 
    
    selected_protein_list <- union(union(y$`AUC > 0.7`, y$`NPX diff > 1`), y$`adj pval < 1E-6`)
    
    # Filter and calculate correlation
    medeca_train_filtered <- 
      medeca_train |> 
      select(Sample, Cancer, selected_protein_list)
    
    cor_selected_proteins <- 
      medeca_train_filtered |> 
      select(-Cancer) |> 
      column_to_rownames("Sample") |> 
      cor(method = "spearman", use = "complete.obs")
    
    cor_plot <- 
      cor_selected_proteins |> 
      pheatmap()
    
    medeca_test_filtered <- 
      medeca_test |> 
      select(Sample, Cancer, selected_protein_list)
    
    # # Run lasso pipeline with selected proteins
    # ml_medeca_filter<- 
    #   do_lasso(variable = "Cancer",
    #            case = "Yes", 
    #            split_train = medeca_train_filtered, 
    #            split_test = medeca_test_filtered) 
    # 
    # #saveRDS(ml_medeca_filter, savepath_data("ML", "ml_medeca_filter.rds"))
    # 
    # protein_importance <- 
    #   ml_medeca_filter$important_proteins |> 
    #   filter(Importance > 0)
    # 
    
    
    
    return(list("selected_proteins_euler" = euler_plot,
                "selected_proteins" = selected_protein_list,
                "selected_proteins_correlation" = cor_plot))#,
    #"protein_importance" = protein_importance,
    # "performance" = ml_medeca_filter$performance)
    
  }

# Function to join data from HPA and MEDECA
combine_data <- 
  function(name_alvez, name) {
    
    data <- 
      de_hpa |> 
      filter(Disease == name_alvez) |> 
      select(Assay, logFC_Alvez = logFC, adj_pval_Alvez = adj.P.Val) |> 
      left_join(de_cancers |> 
                  filter(Cancer == name) |> 
                  select(Assay, logFC, adj_pval = adj.P.Val), 
                by = "Assay") |> 
      mutate(p_value = case_when(
        adj_pval < 0.05 & adj_pval_Alvez < 0.05 ~ "Both studies",
        adj_pval < 0.05 ~ "One study",
        adj_pval_Alvez < 0.05 ~ "One study",
        TRUE ~ "None")) |> 
      mutate(Cancer = name)
    
    return(data)
    
  }

