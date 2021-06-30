library(RColorBrewer)
library(ggthemes)

Axis_themes <- theme(plot.title = element_text(size = 8),
                     axis.title = element_text(size = 7), 
                     axis.text = element_text(size = 6),
                     axis.text.x = element_text(size = 6),
                     legend.text = element_text(size =6),
                     legend.title = element_text(size = 8),
                     strip.text.x = element_text(size = 8))

remove_grid = theme_bw() + theme( panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank())



UMAP_theme <- theme_bw()+theme(axis.text.y = element_blank(), 
                               axis.text.x = element_blank(), 
                               axis.ticks.x= element_blank(),
                               axis.ticks.y= element_blank(),
                               axis.title.x= element_blank(),
                               axis.title.y= element_blank(),
                               strip.text.x = element_text(size = 8),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.border = element_blank(),
                               panel.background = element_rect(colour = "black", fill = NA,size=.6),
                               legend.position = "none",
                               plot.title = element_text(hjust = 0.5, size = 7))


tissue_palette = brewer.pal(12,'Paired')[c(2, 1, 6, 5)]
names(tissue_palette) = c('Diseased_D', 'Remission_D', 'Diseased_E', 'Remission_E')