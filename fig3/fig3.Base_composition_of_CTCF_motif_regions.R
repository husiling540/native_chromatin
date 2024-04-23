library(ggplot2)
library(stringr)
library(RColorBrewer)
library(ComplexHeatmap)
library(pheatmap)

base <- read.table(file = "fig3.Base_composition_of_CTCF_motif_regions.txt", header = FALSE)
colnames(base) <- c("motif_chr", "motif_start", "motif_end", "motif_index", "salt_tolerance", "strand", "peak_group",
                    "Usequence", "Csequence", "Dsequence", "U", "C", "D", "xchip_signal")

base <- base[order(base$salt_tolerance, base$xchip_signal, decreasing = TRUE),]

base_salt_tolerance <- matrix(base$salt_tolerance)

base_matrix <- data.frame(sequence = paste0(base$Usequence, "NNN", base$Csequence, "NNN", base$Dsequence))

base_matrix[, paste0("base", 1:44)] <- str_split_fixed(base_matrix$sequence, "", 44)

rownames(base) <- paste(base$peak_group, base$motif_index, sep = "_")

base_matrix <- base_matrix[, paste0("base", 1:44)]
base_matrix[base_matrix == "A"] <- 1
base_matrix[base_matrix == "T"] <- 3
base_matrix[base_matrix == "C"] <- 5
base_matrix[base_matrix == "G"] <- 7
base_matrix[base_matrix == "N"] <- 9

base_matrix <- data.frame(lapply(base_matrix, as.numeric))
rownames(base_matrix) <- rownames(base)

annotation_row <- data.frame(
  D = factor(base$D),
  U = factor(base$U),
  peak_cluster = factor(base$peak_group),
  salt = base_salt_tolerance
)
rownames(annotation_row) <- rownames(base_matrix)

breaks <- c(0, 2, 4, 6, 8, 10)
colors <- c("#00FB00", "#FF0000", "#0039FF", "#FFFA00", "white")
salt_color = c(colorRampPalette(c("#F6EBE7", "#CF3F23"))(100), colorRampPalette(c("#CF3F23", "#a4161a"))(3000))
row_colors <- list(U = c(`0` = "#F2F2F2", `1` = "red"),
                   D = c(`0` = "#F2F2F2", `1` = "red"),
                   peak_cluster = c(N1 = "#a8e0a8", N2 = "#9f79d3", N3 = "#386CAF", XChIP = "#fc9533"),
                   salt = salt_color)

pdf(file = "fig3.Base_composition_of_CTCF_motif_regions.pdf", width = 6, height = 5.3)
pheatmap(base_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         breaks = breaks,
         color = colors,
         legend = TRUE,
         border_color = NA,
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_colors = row_colors,
         annotation_row = annotation_row)
dev.off()

base_salt_tolerance <- matrix(base$salt_tolerance)
pdf(file = "fig3.Base_heatmap_S-score_legend.pdf", width = 1, height = 5.3)
ComplexHeatmap::pheatmap(base_salt_tolerance,
                         cluster_rows = FALSE,
                         cluster_cols = FALSE,
                         color = salt_color,
                         legend = T,
                         legend_breaks = c(-0.3, 1, 5, 25),
                         border_color = NA,
                         show_rownames = FALSE,
                         show_colnames = FALSE)
dev.off()