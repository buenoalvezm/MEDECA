#### Title: Functions for data analyses 
#### Author: María Bueno Álvez
#### Description: script collecting functions to perform data analyses (DE & ML)
#### Last edited : 12/08/2024

# Data analyses packages
library(limma)
library(tidymodels)
library(themis)

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
  proteins <- 
    data |> 
    select(-Sample, -Cancer) |> 
    colnames()
  
  auc_res <- 
    map_df(proteins, function(protein) {
    protein_exp <- 
      data |> 
      pull(sym(protein))
    
    cancer_class <- 
      data |> 
      pull(Cancer)
    
    res <- pROC::roc(cancer_class, protein_exp)
    
    tibble(Protein = protein,
           AUC = round(as.numeric(res$auc), 2))
    
  })
  
  plot_auc <- 
    auc_res |> 
    ggplot(aes(AUC)) +
    geom_histogram() +
    theme_hpa()
  
  auc_res_cutoff <- 
    auc_res |> 
    filter(AUC > cutoff)
  
  return(list(full_list = auc_res,
              auc_distribution = plot_auc,
              final_list = auc_res_cutoff))
  
  }


# Function to perform lasso analyses
do_lasso <-  
  function(variable,
           case, 
           split_train, 
           split_test#,
           #auc_filter = F
           ) {
    
    # Prepare data - make custom split for current variable
    training_dat <- 
      split_train |> 
      mutate(Class = case_when(!!sym(variable) == case ~ "Case",
                               !!sym(variable) != case ~ "Control")) |> 
      mutate(Class = factor(Class)) |> 
      select(-!!sym(variable)) |> 
      mutate(Sample = as.character(Sample))
    
    testing_dat <- 
      split_test  |> 
      mutate(Class = case_when(!!sym(variable) == case ~ "Case",
                               !!sym(variable) != case ~ "Control")) |> 
      mutate(Class = factor(Class)) |> 
      select(-!!sym(variable))|> 
      mutate(Sample = as.character(Sample))
    
    
    ml_split_custom <- make_splits(training_dat, testing_dat)
    
    
    # # Recipe with ML steps
    # if(auc_filter == T) {
    #   
    #   # Find proteins to retain
    #   proteins <- 
    #     step_filter_auc(data = training_dat, 
    #                     cutoff = 0.75)
    #   
    #   # Keep proteins in training and testing data
    #   training_dat <- 
    #     training_dat |> 
    #     select(Sample, Class, proteins)
    #   
    #   testing_dat <- 
    #     testing_dat |> 
    #     select(Sample, Class, proteins)
    #   
    #   
    # } else {
    #   next
    # }
    # 
    
    # Define recipe
    ml_recipe <- 
      recipe(Class ~ ., data = training_dat) |> 
      step_relevel(Class, ref_level = "Control") |> 
      update_role(Sample, new_role = "id") |> 
      step_normalize(all_numeric_predictors()) |> 
      step_nzv(all_numeric_predictors()) |> 
      #step_corr(all_numeric_predictors()) |> 
      step_impute_knn(all_numeric_predictors()) 
    
    # LASSO model specifications
    glmnet_specs <- 
      logistic_reg() |> 
      set_mode("classification") |> 
      set_engine("glmnet") |> 
      set_args(penalty = tune(), 
               mixture = tune()) #1 
    
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
      grid_latin_hypercube(size = 30)
    
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
      select_best(glmnet_res, metric = "roc_auc", event_level = "second") |> 
      select(-.config)
    
    #Finalize the workflow and fit the final model
    glmnet_wflow <- 
      glmnet_wflow |>  
      finalize_workflow(best_glmnet)
    
    final_glmnet_fit <- last_fit(glmnet_wflow,
                                 ml_split_custom, 
                                 #event_level = "second",
                                 metrics = eval_metrics) 
    
    # Extract model performance
    performance <- 
      final_glmnet_fit |> 
      collect_metrics(event_level = "second") |> 
      select(-.config, -.estimator)
    
    glmnet_auc <- 
      final_glmnet_fit |> 
      collect_metrics(event_level = "second") |> 
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
                "final_workflow" = glmnet_wflow,
                "final_fit" = final_glmnet_fit,
                "predictions" = predictions,
                "performance" = performance,
                "confusion_matrix" = cm,
                "roc_curve" = roc, 
                "important_proteins" = important_proteins))
  }


# Function to join data from HPA and MEDECA
combine_data <- 
  function(name_alvez, name) {
    
    data <- 
      de_pancancer |> 
      filter(Cancer == name_alvez) |> 
      select(Assay, logFC_Alvez = NPX_difference, adj_pval_Alvez = p.adjusted) |> 
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

