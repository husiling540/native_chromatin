library(ggplot2)

native <- read.table(file = "Ratio_of_native_sites.txt", sep = "\t", header = FALSE)
nonnative <- read.table(file = "Ratio_of_nonnative_sites.txt", sep = "\t", header = FALSE)

ratio_stat <- rbind(native, nonnative)
colnames(ratio_stat) <- c("cell_type", "peak_type", "ratio")
ratio_stat$ratio <- round(as.numeric(ratio_stat$ratio) * 100, 1)

ratio_stat$peak_type <- factor(ratio_stat$peak_type, levels = c("native", "non-native"))
ratio_stat$cell_type <- factor(ratio_stat$cell_type, levels = c("0-10cell", "10-20cell", "20-30cell", "30-40cell", "40-50cell", "50-60cell", "60-70cell"))

ggplot(ratio_stat, aes(x = peak_type, y = ratio, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.9, position = "stack") +
  theme_bw() +
  scale_fill_manual(name = "Cell types", values = c("#5A8BBA", "#73B9D1", "#C9DBF1", "#DBEDF3", "#F5DD7A", "#EC8C49", "#BE372D"),
                    labels = c("0-10", "11-20", "21-30", "31-40", "41-50", "51-60", "61-66")) +
  labs(x = NULL, y = "% CTCF sites", title = NULL) +
  geom_text(aes(label = paste0(ratio)), position = position_stack(vjust = 0.5), size = 5, colour = "black") +
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
        plot.title = element_text(colour = "black", face = "bold", size = 16, vjust = 2, hjust = 0)) +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 25, 50, 75, 100), labels = c("0", "25", "50", "75", "100")) +
  theme(legend.position = "right", legend.title = element_text(size = 15), legend.text = element_text(size = 15))

ggsave("fig5.Percentage_of_CTCF_sites_shared_by_multiple_cell_types.pdf", width = 4.2, height = 3.3)