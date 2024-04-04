
# Run differential expression using limma
de_limma_disease <-
  function(data_wide, 
           metadata,
           disease,
           correct = T) {
    
    # Filter data (cohort, variable is 0 or 1)
    dat <-
      data_wide %>% 
      inner_join(metadata %>% 
                   select(Sample, Sex, Age, BMI, Disease), by = "Sample") %>% 
      rename(Group = Disease) %>% 
      mutate(Group = ifelse(Group == disease, "1_Case", "0_Control")) 
    
    # Design a model - add Group, and Sex, Age, BMI
    if(correct == T) {
      design <- model.matrix(~0 + as.factor(dat$Group) + as.factor(dat$Sex) + dat$Age ) #+ dat$BMI
      colnames(design) <- c("control", "case",  "Sex", "Age") #, "BMI"
    } else {
      design <- model.matrix(~0 + as.factor(dat$Group))
      colnames(design) <- c("control", "case")
    }
    
    # Make contrast
    contrast <- makeContrasts(Diff = case - control, levels = design)
    
    dat_fit <- 
      dat %>% 
      select(-Sex, -Age, -BMI, -Group)  %>% 
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
      mutate(Disease = disease)
    
    return(DE_res)
  }
