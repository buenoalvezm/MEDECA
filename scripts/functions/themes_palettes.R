#### Title: Themes & palettes
#### Author: María Bueno Álvez
#### Description: script collecting themes & palettes
#### Last edited : 12/08/2024

# Packages for themes & palettes
library(RColorBrewer)
library(paletteer)

# Levels 
cancer_groups <- 
  c("Lymfom",
    "Myelom", 
    "Lung",
    "Neuroendocrine tumors",
    "Pancreas or gall bladder or bile duct",
    "Breast",
    "Ovarian or tubar",
    "No diagnosis")

# Palettes
pal_panel <- c("#28C7E3",
               "#FA1E00",
               "#0A56A1",
               "#FDC600")
pal_de <-
  c("not significant" = "#D3D3D3",
    "significant up" = "#FF7176",
    "significant down" = "#92C9DA")

pal_cancer <- c("Yes" = "#E16C54", "No" = "#A0CBC7")

pal_sex <- c("Male" = "#3E4A57", "Female" = "#CEC6B3")

pal_disease_type<- 
  c("Cancer" = "#E16C54", 
    "No diagnosis" = "#BBDABB", 
    "Autoimmune" = "#E9E7AB", 
    "Infectious" = "#E9E7AB", 
    "Inflammatory" = "#E9E7AB", 
    "Other" = "#E9E7AB")

pal_category <- 
  c("Cancer" = "#E16C54", 
    "No diagnosis" = "#BBDABB",
    "Other diseases" = "#F4B183")


pal_cancers <- c("#08585A",
                 "#66C2A5",
                 "#ADC74F",
                 "#FDB36A",
                 "#D7485A",
                 "#E8A29A",  
                 "#603479",
                 "grey")
names(pal_cancers) <- cancer_groups

pal_cancers_overlap <- c("Lymfom" = "#08585A",
                         "Myelom" = "#66C2A5",
                         "Lung" = "#ADC74F",
                         "Breast" = "#E8A29A",  
                         "Ovarian tubar" = "#603479")

pal_red <- rev(colorRampPalette(brewer.pal(6, "Reds"))(10))

pal_binary <- c("Yes" = "red",
                "No" = "grey")

pal_heatmap <-
  brewer.pal(9, name = "YlOrRd") |>  
  colorRampPalette()

pal_controls <- 
  c("No diagnosis" = "#BBDABB",
    "Infectious" = "#E5B9AD", 
    "Autoimmune" = "#D98994",
    "Inflammatory" = "#D0587E",
    "Other" = "grey")

pal_pvalue <- c("Both studies" = "#945785", 
                "One study" = "#DAA4AB", 
                "None" = "grey")

# HPA theme
theme_hpa <- 
  function(angled = F, axis_x = T, axis_y = T, facet_title = T) {
    t <- 
      theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border = element_blank(),
        plot.title = element_text(face = "bold",
                                  size = rel(1), hjust = 0.5),
        plot.subtitle=element_text(face = "bold",hjust = 0.5, size=rel(1),vjust=1),
        axis.title = element_text(face = "bold",size = rel(1)),
        axis.ticks.length = unit(.25, "cm"),
        axis.line = element_line(linewidth = 0.5),
        axis.text = element_text(size = rel(1), color = 'black'),
        legend.key = element_blank(),
        legend.position = "right",
        legend.text = element_text(size=rel(0.8)),
        legend.key.size= unit(0.7, "cm"),
        legend.title = element_text(size=rel(1)),
        plot.margin=unit(c(10,5,5,5),"mm"),
        strip.background=element_rect(colour="grey90",fill="grey90"),
        strip.text = element_text(face="bold")
      )
    
    if(angled) {
      t <- t + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
    }
    
    if(axis_x == F) {
      t <- t +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              axis.title.x = element_blank())
    } 
    
    if(axis_y == F) {
      t <- t +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.line.y = element_blank(),
              axis.title.y = element_blank())
    }
    if(facet_title == F) {
      t <- t + theme(strip.text = element_blank())
    }
    return(t)
  }



