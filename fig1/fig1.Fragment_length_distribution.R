library(dplyr)
library(reshape2)
library(ggplot2)

fragment_length <- data.frame(fragment_length = c(10:250))

for (i in list.files(pattern = ".txt")) {
  short_name <- sub(".txt", "", sub("readlength_", "", i))
  expeiment <- read.table(file = i, header = FALSE)
  expeiment_nrow <- nrow(expeiment)
  expeiment_t <- data.frame(table(expeiment))
  colnames(expeiment_t) <- c("fragment_length", short_name)
  expeiment_t[, short_name] <- expeiment_t[, short_name] / expeiment_nrow
  fragment_length$fragment_length <- as.character(fragment_length$fragment_length)
  fragment_length <- fragment_length %>%
    left_join(expeiment_t, by = "fragment_length")
}

fragment_length[is.na(fragment_length)] <- 0

fragment_length <- fragment_length[, c("fragment_length", "loMNase", "DNase", "ATAC", "MNase_304U", "MNase_79.2U", "MNase_20.6U", "MNase_5.4U")]
fragment_length$fragment_length <- as.character(fragment_length$fragment_length)
fragment_length_m <- melt(fragment_length, ID = c("fragment_length"))
fragment_length_m$fragment_length <- as.integer(fragment_length_m$fragment_length)

ggplot(fragment_length_m, aes(x = fragment_length, y = variable)) +
  geom_tile(aes(fill = value)) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_bw() +
  scale_fill_gradient("Fraction", low = "white", high = "#32469B") +
  scale_x_continuous(expand = c(0, 0), labels = c(10, "", "", "", 50, rep("", 4), 100, rep("", 4), 150, rep("", 4), 200, rep("", 4), 250, rep("", 4), 300), breaks = seq(10, 300, 10)) +
  scale_y_discrete(position = "left", labels = c("loMNase", "DNase", "ATAC", "MNase 304U", "MNase 79.2U", "MNase 20.6U", "MNase 5.4U")) +
  theme(legend.position = "right",
        legend.box = "horizontal",
        legend.key.size = unit(15, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
        legend.text = element_text(colour = "black", size = 15),
        panel.border = element_blank(),
        axis.title.x = element_text(colour = "black", size = 15, hjust = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 14, colour = "black"),
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 14, colour = "black"),
        axis.line = element_blank(),
        axis.ticks.x = element_line(size = 0.5, color = "black"),
        axis.ticks.y = element_blank(),
        text = element_text(size = 15, family = "sans", colour = "black"))

ggsave('fig1.Fragment_length_distribution.pdf', width = 5.5, height = 2.5, dpi = 600)