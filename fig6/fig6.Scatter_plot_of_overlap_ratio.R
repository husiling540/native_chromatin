library(ggplot2)
library(dplyr)

ratio_loss <- read.table(file = "Overlap_ratio_of_K562_loMNase_lost_sites.sim.txt", sep = "\t", header = FALSE)
colnames(ratio_loss) <- c("TF", "peak_count", "lost_ratio")

ratio_stable <- read.table(file = "Overlap_ratio_of_K562_loMNase_stable_sites.sim.txt", sep = "\t", header = FALSE)
colnames(ratio_stable) <- c("TF", "peak_count", "stable_ratio")

overlap_ratio <- ratio_loss %>% left_join(ratio_stable, by = "TF")

triangle <- data.frame(x = c(0, 0, 1), y = c(0, 1, 1))

ggplot(overlap_ratio) +
  geom_point(aes(x = stable_ratio, y = lost_ratio, color = TF), alpha = 1, size = 3) +
  geom_polygon(data = triangle, aes(x = x, y = y), fill = "grey", alpha = 0.2) +
  labs(title = "", x = 'Overlap ratio of stable sites', y = "Overlap ratio of lost sites") +
  theme_bw() +
  theme(legend.box = "horizontal",
        legend.key.size = unit(18, "pt"),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
        legend.text = element_text(colour = "black", size = 13),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        axis.title.x = element_text(colour = "black", size = 14, vjust = 0),
        axis.title.y = element_text(colour = "black", size = 14, vjust = 2),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 13, colour = "black"),
        axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 13, colour = "black")) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual("", values = c("grey", "#DB453E", "black", "#386CAF", "#8DC73F", "#EE9B41", "#AC5895", "#81CFD9", "#8D81AF", "#C29B39")) +
  scale_shape(guide = guide_legend(title.position = "top", size = 13)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(-1, 1, 0.2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(-1, 1, 0.2), expand = c(0, 0)) +
  coord_fixed(ratio = 1)

ggsave("fig6.Scatter_plot_of_overlap_ratio.pdf", width = 4.5, height = 3.4)