library(ggplot2)
library(ggpubr)

nonSUMO_Sscore <- read.table(file = "K562_CTCF_non-SUMO_sites.salt_tolerance_score.txt", header = F)
colnames(nonSUMO_Sscore) <- c("site_chr", "site_start", "site_end", "motif_chr", "motif_start", "motif_end", "motif_index", "salt_tolerance", "strand", "75mM", "150mM", "225mM")
nonSUMO_Sscore$peak_cluster <- "nonSUMO"

SUMO_Sscore <- read.table(file = "K562_CTCF_SUMO_sites.salt_tolerance_score.txt", header = F)
colnames(SUMO_Sscore) <- c("site_chr", "site_start", "site_end", "motif_chr", "motif_start", "motif_end", "motif_index", "salt_tolerance", "strand", "75mM", "150mM", "225mM")
SUMO_Sscore$peak_cluster <- "SUMO"

motif_Sscore <- rbind(SUMO_Sscore, nonSUMO_Sscore)

motif_Sscore$peak_cluster <- factor(motif_Sscore$peak_cluster, levels = c("SUMO", "nonSUMO"))

my_comparisons <- list(c("SUMO", "nonSUMO"))

ggplot(motif_Sscore, aes(x = peak_cluster, y = salt_tolerance, fill = peak_cluster)) +
  geom_boxplot(width = 0.8, outlier.shape = NA) +
  scale_fill_manual(values = c("#AD8A75", "#96AFD2"), labels = c("SUMO", "nonSUMO")) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.7, vjust = 1, colour = "black", size = 18, angle = 30),
        axis.text.y = element_text(size = 18, face = "plain", colour = "black"),
        axis.title.x = element_text(size = 18, face = "plain", colour = "black", vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 18, face = "plain", colour = "black", vjust = 1, hjust = 0.5),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks = element_line(colour = "black", size = 0.5),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
        plot.title = element_text(colour = "black", face = "plain", size = 18, vjust = 8, hjust = 1)) +
  coord_cartesian(ylim = c(-1, 9)) +
  scale_x_discrete(breaks = c("SUMO", "nonSUMO"), labels = c("SUMO2/3", "non-SUMO2/3")) +
  scale_y_continuous(breaks = seq(0, 10, 3), labels = seq(0, 10, 3)) +
  labs(title = NULL, x = NULL, y = "S-score") +
  stat_compare_means(comparisons = my_comparisons, label.y = c(7.5), method = "wilcox.test", label = "p.signif")

ggsave("fig4.Boxplot_of_S-score.pdf", width = 3, height = 4)