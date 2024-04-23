library(reshape2)
library(ggplot2)

motif_percent <- read.table(file = "Number_of_peaks_containing_CTCF_motif.txt", sep = "\t", header = FALSE)
colnames(motif_percent) <- c("exp", "split", "_with_motif", "_without_motif", "ratio")

motif_percent_m <- melt(motif_percent, id.vars = c("exp", "split", "ratio"))
colnames(motif_percent_m) <- c("exp", "split", "ratio", "peak_group", "count")

motif_percent_m$peak_group <- factor(motif_percent_m$peak_group, levels = c("_without_motif", "_with_motif"))

motif_percent_m$exp_split <- paste0(motif_percent_m$exp, "_", motif_percent_m$split)
motif_percent_m$exp <- factor(motif_percent_m$exp, levels = c("NChIP-public", "CUT-Tag", "CUT-RUN", "NChIP"))

motif_percent_m$split <- factor(motif_percent_m$split, levels = c(1:3), labels = c("0.7 M", "1.4 M", "2.1 M"))

motif_percent_m$exp_peak <- paste0(motif_percent_m$exp, motif_percent_m$peak_group)
motif_percent_m$exp_peak <- factor(motif_percent_m$exp_peak, levels = c("NChIP-public_without_motif", "CUT-Tag_without_motif", "CUT-RUN_without_motif", "NChIP_without_motif", "NChIP-public_with_motif", "CUT-Tag_with_motif", "CUT-RUN_with_motif", "NChIP_with_motif"))

ggplot(data = motif_percent_m, aes(x = exp, y = count, fill = exp_peak, color = exp_peak)) +
  geom_bar(stat = "identity", width = 0.9, size = 1, position = "stack") +
  scale_fill_manual(name = NULL, values = c("white", "white", "white", "white", "#7FC87F", "#C29B39", "#386CAF", "#A67BB0")) +
  scale_color_manual(name = NULL, values = c("#7FC87F", "#C29B39", "#386CAF", "#A67BB0", "#7FC87F", "#C29B39", "#386CAF", "#A67BB0")) +
  labs(x = NULL, y = "Number of peaks", title = NULL) +
  facet_wrap(~split) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        plot.margin = unit(c(0.5, 0.2, 0.2, 0.2), "lines"),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent"),
        axis.ticks = element_line(color = "black", size = 0.5),
        axis.ticks.x = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 28, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 28, colour = "black", angle = 90, hjust = 0.5, vjust = 0.5),
        axis.title.x = element_text(size = 28, color = "black", hjust = 0.5),
        axis.title.y = element_text(size = 28, color = "black", hjust = 0.5, angle = 90),
        plot.title = element_text(colour = "black", face = "bold", size = 28, vjust = 2, hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(size = 28, angle = 0),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 28)) +
  scale_y_continuous(position = "left", limits = c(0, 20000), expand = c(0, 0), breaks = c(0, 5e3, 10e3, 15e3, 20e3), labels = format(c(0, 5e3, 10e3, 15e3, 20e3), big.mark = ",")) +
  scale_x_discrete(labels = c("NChIP", "CUT-RUN", "CUT-Tag", "NChIP"))

ggsave('fig2.Numbers_of_peaks_containing_CTCF_motif.pdf', width = 12, height = 9.5, dpi = 600)