library(ggplot2)
library(ggpubr)

load_data <- function(filename, sample_name) {
  data <- read.table(filename, header = FALSE, sep = "\t")
  names(data)[1:4] <- c("chr", "start", "end", "length")
  data$sample <- sample_name
  return(data)
}

loMNase <- load_data("violin_K562_loMNase_fragment_length.txt", "loMNase")
DNase <- load_data("violin_K562_DNase_fragment_length.txt", "DNase")
ATAC <- load_data("violin_K562_ATAC_fragment_length.txt", "ATAC")

All_width <- rbind(loMNase, DNase, ATAC)

ggplot(All_width, aes(x = sample, y = length)) +
  geom_violin(trim = FALSE, color = NA, aes(fill = sample), size = 0.1, width = 0.9, scale = "area", alpha = 0.8) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), size = 0.5, alpha = 1, outlier.shape = NA, fill = NA) +
  scale_fill_manual(values = c("#2A5CA3", "#C36518", "grey")) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, colour = "black", size = 17, angle = 20),
        axis.text.y = element_text(size = 18, face = "plain", colour = "black"),
        axis.title.x = element_text(size = 18, face = "plain", colour = "black", vjust = -1.5, hjust = 0.5),
        axis.title.y = element_text(size = 18, face = "plain", colour = "black", vjust = 2, hjust = 0.5),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.5), "lines"),
        plot.title = element_text(colour = "black", face = "plain", size = 18, vjust = 8, hjust = 1)) +
  labs(title = "", x = "", y = "Peak-containing\nfragment size (bp)", size = 19) +
  scale_x_discrete(breaks = c("loMNase", "DNase", "ATAC"), labels = c("loMNase", "DNase", "ATAC")) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 600, 100), labels = seq(0, 600, 100)) +
  coord_cartesian(ylim = c(0, 520)) +
  stat_compare_means(comparisons = list(c("loMNase", "DNase"), c("loMNase", "ATAC")), label.y = c(200, 450), method = "wilcox.test", label = "p.signif", paired = FALSE, alternative = "two.sided")

ggsave("fig1.Comparing_the_size_difference_of_fragments_containing_peak_summits.pdf", width = 3.8, height = 4.2)