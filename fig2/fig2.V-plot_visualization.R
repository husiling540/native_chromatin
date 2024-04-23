args <- commandArgs(TRUE)

library(ggplot2)
library(scales)

native_v <- read.table(file = paste0(args[1], ".midpoint_to_motif.bed"), sep = "\t", header = F)
colnames(native_v) <- c("read_chr", "read_midstart", "read_midend", "read_length",
                        "motif_chr", "motif_start", "motif_end", "motif_index", "motif_na", "motif_strand", "distance")
native_v <- native_v[(native_v$read_chr != "."), c("read_length", "distance")]
native_v$read_length <- as.numeric(native_v$read_length)

ggplot(native_v, aes(x = distance, y = read_length)) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE, n = 500) +
  scale_fill_gradient2(low = "white", high = "dodgerblue4") +
  theme_bw() +
  theme_classic() +
  theme(panel.border = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.x = element_line(color = "black", size = 0.6),
        axis.ticks.y = element_line(color = "black", size = 0.6),
        axis.title.x = element_text(size = 21, vjust = 0),
        axis.title.y = element_text(size = 21),
        axis.text.x = element_text(size = 21, colour = 'black', angle = 30, hjust = 0.6, vjust = 0.7),
        axis.text.y = element_text(size = 21, colour = 'black'),
        plot.title = element_text(colour = "black", size = 21, hjust = 0.5),
        plot.margin = unit(c(.1, .3, 0.1, .1), 'inches'),
        legend.position = "none") +
  labs(title = args[1], x = "Distance from CTCF motif (bp)", y = "Fragment length (bp)", size = 18) +
  scale_x_continuous(limits = c(-100, 100), expand = c(0, 0), breaks = seq(-100, 100, 50)) +
  scale_y_continuous(limits = c(0, 200), expand = c(0, 0), breaks = seq(0, 200, 50))

ggsave(paste0("fig2_", args[1], "_vplot.pdf"), width = 4, height = 3.5)