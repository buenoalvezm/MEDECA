
# Functions for visualization

## Generate volcano plot from differential expression results                                                                                                                                  
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
