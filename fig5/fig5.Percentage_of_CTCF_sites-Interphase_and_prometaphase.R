library(ggplot2)

mitosis <- read.table(file = "U2OS_mitosis_peak_fraction.txt", sep = "\t", header = F)
interphase <- read.table(file = "U2OS_interphase_peak_fraction.txt", sep = "\t", header = F)

all <- rbind(mitosis, interphase)
colnames(all) <- c("phase", "native_or_non", "phase_peak", "native_nonnative_peak", "overlap_phase_peak", "overlap_ratio")

all$native_or_non <- factor(all$native_or_non, levels = c("nonnative", "native"))

ggplot(all, aes(x = phase, y = overlap_ratio, fill = native_or_non)) +
  geom_bar(stat = "identity", width = 0.9, position = "stack") +
  theme_bw() +
  scale_fill_manual(name = "", values = c("#AEB7C8", "#3F73B2"), labels = c("non-native", "native")) +
  labs(x = NULL, y = "% CTCF sites", title = NULL) +
  geom_text(aes(label = paste0(round(overlap_ratio * 100, 1))), position = position_stack(vjust = 0.5), size = 5.5, colour = "black") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_line(color = "black", size = 0.5),
        axis.line = element_line(color = "black", size = 0.5),
        axis.text.x = element_text(size = 16, colour = "black", angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title.x = element_text(size = 16, color = "black", hjust = 0.5),
        axis.title.y = element_text(size = 16, color = "black", hjust = 0.5),
        plot.title = element_text(colour = "black", face = "bold", size = 16, vjust = 2, hjust = 0),
        legend.position = "top",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)) +
  scale_x_discrete(labels = c("Interphase", "Prometaphase")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = seq(0, 100, 25))

ggsave("fig5.Percentage_of_CTCF_sites-Interphase_and_prometaphase.pdf", width = 3.7, height = 3.5, dpi = 600)