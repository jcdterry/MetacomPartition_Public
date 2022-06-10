species_plot <- function(data, plotMain = NULL, colorVar = NULL, colorLegend = "none"){
  data %>% 
    ggtern(aes(x = env, z = spa, y = codist, size = r2)) +
    geom_point(aes_string(color = colorVar), alpha = 0.8) +
    scale_T_continuous(limits=c(0,1.0),
                       breaks=seq(0,1,by=0.1),
                       labels=seq(0,1,by=0.1)) +
    scale_L_continuous(limits=c(0.0,1),
                       breaks=seq(0,1,by=0.1),
                       labels=seq(0,1,by=0.1)) +
    scale_R_continuous(limits=c(0.0,1.0),
                       breaks=seq(0,1,by=0.1),
                       labels=seq(0,1,by=0.1)) +
    labs(title = plotMain,
         x = "Environment Max",
         xarrow = "Environment",
         y = "Co-distribution Max",
         yarrow = "Co-Distribution",
         z = "Space Max", 
         zarrow = "Spatial Autocorrelation") +
    theme_light() +
    theme_showarrows() +
    #scale_colour_brewer(palette = "Set1") +
    #scale_colour_brewer(palette = "Spectral") +
    #scale_color_viridis_d() +
    scale_size_area(limits = c(0,1), breaks = seq(0,1,0.2)) +
    guides(color = guide_legend(colorLegend, order = 2), 
           size = guide_legend(title = expression(R^2), order = 1)) +
    theme(panel.grid = element_line(color = "darkgrey"),
          axis.title = element_text(size = 8))
}
