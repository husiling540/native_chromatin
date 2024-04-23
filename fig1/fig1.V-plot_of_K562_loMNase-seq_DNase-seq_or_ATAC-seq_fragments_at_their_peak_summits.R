args <- commandArgs(TRUE)

library(ggplot2)
library(scales)

v_data <- read.table(file = paste0("sample.K562_", args[1], ".midpoint_to_", args[1], "_summits.bed"), sep = "\t", header = FALSE)
colnames(v_data) <- c("read_chr", "read_midstart", "read_midend", "read_length", "motif_chr", "motif_start", "motif_end", "motif_index", "motif_na", "motif_strand", "distance")
v_data <- v_data[v_data$read_chr != ".", c("read_length", "distance")]
v_data$read_length <- as.numeric(v_data$read_length)

v_data_sub1 <- v_data[(v_data$read_length - abs(2 * v_data$distance) >= 0) & (v_data$read_length <= 70),]
v_data_sub2 <- v_data[(abs(v_data$distance) <= 100) & (v_data$read_length <= 100),]

ratio <- paste0(round(nrow(v_data_sub1) / nrow(v_data_sub2) * 100, 0), "%")

ggplot(v_data, aes(x = distance, y = read_length)) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE, n = 500) +
  geom_abline(intercept = 0, slope = -2, color = "#C2212B") +
  geom_abline(intercept = 0, slope = 2, color = "#C2212B") +
  geom_hline(yintercept = c(70), color = "#C2212B") +
  scale_fill_gradient2(low = "white", high = "dodgerblue4") +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black", size = 0.6),
        axis.title = element_text(size = 21),
        axis.text = element_text(size = 21, colour = 'black'),
        plot.title = element_text(colour = "black", size = 21, hjust = 0.5),
        plot.margin = unit(c(0.1, 0.3, 0.1, 0.1), 'inches'),
        legend.position = "none") +
  labs(title = args[1], x = "Distance from peak summit (bp)", y = "Fragment length (bp)", size = 18) +
  scale_x_continuous(limits = c(-100, 100), expand = c(0, 0), breaks = seq(-100, 100, 50)) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0), breaks = seq(0, 100, 50)) +
  annotate('text', x = 28, y = 10, label = ratio, size = 6, colour = "#C2212B")

ggsave(paste0("fig1.K562_", args[1], ".midpoint_to_", args[1], "_summits.pdf"), width = 4, height = 3.5)