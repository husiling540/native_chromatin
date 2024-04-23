library(ggplot2)
library(ggseqlogo)

m <- read.table(file = "meme_K562_native_sites_6library_noCTCFmotif_MOTIF_1.txt", skip = 3)

colnames(m) <- c("T", "G", "C", "A")

m <- m[seq(nrow(m), 1, by = -1),]

ggseqlogo(t(m), method = 'bits', seq_type = "dna") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2), breaks = c(0, 1, 2)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq_len(nrow(m))) +
  labs(x = NULL) +
  theme(axis.title = element_text(size = 15, colour = 'black'),
        text = element_text(size = 15, family = "sans", colour = "black"),
        axis.ticks.y = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(linetype = 1, color = "black", size = 0.5),
        axis.text.y = element_text(size = 15, colour = 'black'),
        axis.text.x = element_text(size = 15, colour = 'black')) +
  coord_fixed(ratio = 3)

ggsave("fig3.CTCF_V-motif.pdf", width = 6.2, height = 2)