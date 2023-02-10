library(ggplot2)

myaxis <- 
  #axis
  theme(axis.title = element_text(size = 15),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 15, color = "black")) 
myplot <-
  #plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  #legend
  theme(legend.key = element_rect(fill = "white", colour = "white"),
        legend.position = "right", 
        legend.text = element_text( color = "black", size = 15),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, 'cm'))  +
  theme(plot.title = element_text(size = 15, hjust = 0.5, colour = "black"))


