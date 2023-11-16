library(ggplot2)
args <- commandArgs(TRUE)

# Set working directory
setwd("/mnt/disk3/husiling/private/data/native_chromatin_summary/Fig1/vplot")

# Read data from file
vdata <- read.table(file = paste0("K562_loMNase.midpoint_to_", args[1], ".bed"), sep = "\t", header = FALSE)
colnames(vdata) <- c("read_chr", "read_midstart", "read_midend", "read_length",
                        "motif_chr", "motif_start", "motif_end", "motif_index", "motif_na", "motif_strand", "distance")
vdata <- vdata[(vdata$read_chr != "."), c("read_length", "distance")]
vdata$read_length <- as.numeric(vdata$read_length)

#plot
ggplot(vdata, aes(x = distance, y = read_length)) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE, n = 500) +
  scale_fill_gradient2(low = "white", high = "dodgerblue4") +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black", size = 0.6),
        axis.title = element_text(size = 21),
        axis.text.x=element_text(size=21, colour = 'black', angle = 30, hjust = 0.6, vjust = 0.7),
        axis.text.y=element_text(size=21, colour = 'black'),
        plot.title = element_text(colour = "black", size = 21, hjust = 0.5),
        plot.margin = unit(c(.1, .3, 0.1, .1), 'inches'),
        legend.position = "none") +
  labs(title = args[1], x = "Distance from motif center (bp)", y = "Fragment length (bp)", size = 18) +
  scale_x_continuous(limits = c(-100, 100), expand = c(0, 0), breaks = seq(-100, 100, 50)) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0), breaks = seq(0, 100, 50))

# Save plot as PDF
ggsave(paste0("/mnt/disk3/husiling/private/data/native_chromatin_summary/Fig1/vplot/Fig1.4_K562_loMNase.midpoint_to_", args[1], "_vplot_test.pdf"), width = 5, height = 4)