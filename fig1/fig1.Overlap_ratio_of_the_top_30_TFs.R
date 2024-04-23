library(ggplot2)

overlap_ratio <- read.table(file = "loMNase_overlap_ratio.sim.txt", sep = "\t", header = F)
colnames(overlap_ratio) <- c("TF", "experiment", "peak_count", "ratio")

overlap_ratio <- overlap_ratio[overlap_ratio$TF %in% head(unique(overlap_ratio$TF), n = 30),]
overlap_ratio$TF <- factor(overlap_ratio$TF, levels = rev(unique(overlap_ratio$TF)))
overlap_ratio <- overlap_ratio[order(overlap_ratio$ratio, decreasing = TRUE),]
overlap_ratio_uniq <- overlap_ratio[!duplicated(overlap_ratio$TF),]

ggplot(data = overlap_ratio, mapping = aes(x = ratio, y = TF, colour = peak_count)) +
  geom_point(size = 2) +
  scale_colour_gradient("Number of \nfactor peaks", low = 'lightblue', high = '#19315F') +
  labs(title = NULL, x = "Overlap ratio", y = NULL) +
  theme_bw() +
  scale_x_continuous(position = "top") +
  theme(legend.position = "right",
        legend.box = "horizontal",
        legend.key.size = unit(15, "pt"),
        plot.margin = unit(c(0.2, 0.2, 0.5, 0.5), "lines"),
        legend.text = element_text(colour = "black", size = 13),
        panel.border = element_rect(fill = NA, size = 0.5, color = 'black'),
        axis.title.x = element_text(colour = "black", size = 13, hjust = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 14, colour = "black"),
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 14, colour = "black"),
        axis.ticks = element_line(size = 0.5, color = "black"),
        text = element_text(size = 13, family = "sans", colour = "black"))

ggsave("fig1.Overlap_ratio_of_the_top_30_TFs.pdf", width = 5, height = 8)