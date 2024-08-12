#### Title: Themes & palettes
#### Author: María Bueno Álvez
#### Description: script collecting themes & palettes
#### Last edited : 12/08/2024

# Packages for themes & palettes
library(RColorBrewer)

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

pal_group <- 
  c("Cancer" = "#E16C54", 
    "Healthy" = "#BBDABB", 
    "Autoimmune" = "#E4A6B0", 
    "Infectious" = "#E3D0A6", 
    "Inflammatory" = "#7A76A9", 
    "Other" = "grey")

pal_group_3 <- 
  c("Cancer" = "#E16C54", 
    "Healthy" = "#BBDABB",
    "Other diseases" = "#E9E7AB")

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

