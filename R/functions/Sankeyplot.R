# Create Sankey plot for sample processing
library(randomcoloR)
library(ggsankey)
library(dplyr)


load(file = "~/git/ped_CapTCRseq/data/tumor_sampleprocessing.RData")

#long format
df <- make_long(ffpe_samples, Shipped, `Library preparation`, `Successful capture`, `Deep sequencing`, `TCR analysis`, Group, Tumor)

# change order so failed samples apprear last
df$next_node <- factor(df$next_node, levels = c('ALCL','BL','BLL','DLBCL','ERMS','HD',
                                                'Lymphoma','n = 4','n = 17','n = 2','n = 23',
                                                'n = 21','n = 25','NB','OS','PMBCL','Solid'))

df$node <- factor(df$node, levels = c('ALCL','BL','BLL','DLBCL','ERMS','HD',
                                      'Lymphoma','n = 4','n = 17','n = 2','n = 23',
                                      'n = 21','n = 25','NB','OS','PMBCL','Solid'))


set.seed(5250)
myColors <- distinctColorPalette(17)
names(myColors) <- levels(factor(df$node))

#failed samples color grey
myColors["n = 2"] <- "grey"
myColors["n = 4"] <- "grey"

# Tumor types white
myColors[names(myColors) %in% c("ALCL", "BL", "BLL", "DLBCL", "ERMS", "HD", "PMBCL", "NB", "OS")] <- "white"

sankey <- ggplot(df, aes(x = x, 
                         next_x = next_x, 
                         node = node, 
                         next_node = next_node,
                         fill = node, color = node,
                         label = node)) +
  geom_sankey(flow.alpha = 0.5, na.rm = TRUE) +
  theme_sankey(base_size = 16) +
  geom_sankey_label(size = 5, hjust = 0.5,
                    color = 1, #color labels black
                    label.size = 0, # remove borders for labels
                    fill = alpha("white",0), #background transparent
                    na.rm = TRUE) +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20, color = "black"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  scale_fill_manual(values = myColors, na.value="white") + #color NAs white
  scale_color_manual(values = myColors, na.value="white") +  #color NAs white
  labs(title = "Tumor sample processing")

pdf("~/git/ped_CapTCRseq/plots/tumor_sankey.pdf",
    width = 20, height = 10)
sankey
dev.off()




