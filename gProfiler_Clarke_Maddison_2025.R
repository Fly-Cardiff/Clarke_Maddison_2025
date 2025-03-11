library("ggplot2")
library("ggrepel")

##### knockdown experiment 
mainAB <- read.csv("/RNASeq/MARKDUP/Elav_KD_WWOX/newg/gProfiler_dmelanogaster_10-03-2025_18-04-32__intersections_ABMain_KD.csv")
mainAB[11,1] <- "TF"

kdAB <- read.csv("/RNASeq/MARKDUP/Elav_KD_WWOX/newg/gProfiler_dmelanogaster_10-03-2025_18-06-20__intersections_KD_AB.csv")
kdAB[11,1] <- "WP"
kdAB[12,1] <- "TF"

kdlacz <- read.csv("/RNASeq/MARKDUP/Elav_KD_WWOX/newg/gProfiler_dmelanogaster_10-03-2025_18-05-27__intersections_KD_LACZ.csv")
kdlacz[18,1] <- "WP"

##### Overexpression experiment 

mainAB <- read.csv("/RNASeq/MARKDUP/Elav_OE_WWOX/gProfiler_dmelanogaster_10-03-2025_17-59-12__intersections_mainAB_OE.csv")
mainAB[42,1] <- "WP"

OEAB <- read.csv("/RNASeq/MARKDUP/Elav_OE_WWOX/gProfiler_dmelanogaster_10-03-2025_18-02-24__intersections_OE_AB.csv")

OELacZ <- read.csv("/RNASeq/MARKDUP/Elav_OE_WWOX/gProfiler_dmelanogaster_10-03-2025_18-00-18__intersections_OE_LacZ.csv")
OELacZ[12,1] <- "WP"
OELacZ[13,1] <- "TF"


########## PLOTTING ##################
# for knockdown exp (limits=c(1,80),breaks=c(1,10,20,40,80)) and ylim 6
# For OE --> ((limits=c(1,400),breaks=c(1,4,40,400)) and y lim 24

toplot <- kdlacz

p <- ggplot(toplot, aes(y=negative_log10_of_adjusted_p_value, x=term_name, size=intersection_size)) + geom_point() + ylab("-Log10 Padj") + labs(size = "N genes")
p2 <- p + facet_grid(~source, space = 'free_x') + 
  theme_bw() + theme(panel.grid.major.x = element_blank(),  # Remove vertical grid lines
                     panel.grid.major.y = element_line(colour = "grey80",linetype = "dashed"),  # Keep horizontal grid lines at tick points
                     panel.grid.minor = element_blank(), strip.text = element_text(family="Arial", face = "bold", colour = "black", size = 12)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(family="Arial", face = "bold", colour = "black", size = 13), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(),legend.position = "bottom", legend.text = element_text(family = "Arial", color="black", size=13), legend.title = element_text(family = "Arial", color="black", size=13), axis.text= element_text(family="Arial", face = "bold", colour = "black", size = 13)) +
  scale_y_continuous(limits = c(0, 24), breaks = seq(0, 24, 3)) +
  scale_size_continuous(limits=c(1,400),breaks=c(1,4,40,400))


ggsave(plot = p2, units = "mm",
  filename = "/RNASeq/MARKDUP/Elav_OE_WWOX/gprof.png",
  width=200,
  height = 90,
  dpi = 600
)



