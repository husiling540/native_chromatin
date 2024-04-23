library(ChIPpeakAnno)

CTCF_75mMNaCl_merge <- read.table("K562_NChIP_CTCF_75mMNaCl_2rep_peak.bed", header = F, sep = "\t")
names(CTCF_75mMNaCl_merge) <- c("chr", "start", "end", "name")

gr_CTCF_75mMNaCl_merge <- GRanges(CTCF_75mMNaCl_merge, format = "BED", header = FALSE)

CTCF_150mMNaCl_merge <- read.table("K562_NChIP_CTCF_150mMNaCl_2rep_peak.bed", header = F, sep = "\t")
names(CTCF_150mMNaCl_merge) <- c("chr", "start", "end", "name")

gr_CTCF_150mMNaCl_merge <- GRanges(CTCF_150mMNaCl_merge, format = "BED", header = FALSE)

CTCF_225mMNaCl_merge <- read.table("K562_NChIP_CTCF_225mMNaCl_2rep_peak.bed", header = F, sep = "\t")
names(CTCF_225mMNaCl_merge) <- c("chr", "start", "end", "name")

gr_CTCF_225mMNaCl_merge <- GRanges(CTCF_225mMNaCl_merge, format = "BED", header = FALSE)

ol <- findOverlapsOfPeaks(gr_CTCF_75mMNaCl_merge, gr_CTCF_150mMNaCl_merge, gr_CTCF_225mMNaCl_merge)

pdf("fig3.VennDiagram.pdf", width = 9, height = 9)
makeVennDiagram(ol, fill = c("#A085BD", "#82B2E2", "#E0853F"), col = NA, cat.col = 'black', cat.cex = 4, cex = 4, lwd = 0, NameOfPeaks = c("N75mM", "N150mM", "N225mM"), margin = 0.1, fontfamily = "sans")
dev.off()

N1 <- data.frame(ol$peaklist$`gr_CTCF_75mMNaCl_merge///gr_CTCF_150mMNaCl_merge///gr_CTCF_225mMNaCl_merge`)
N2 <- data.frame(ol$peaklist$`gr_CTCF_75mMNaCl_merge///gr_CTCF_150mMNaCl_merge`)
N3 <- data.frame(ol$peaklist$gr_CTCF_75mMNaCl_merge)

write.table(N1, file = "N1_sites.bed", col.names = F, row.names = F, sep = "\t", quote = F)
write.table(N2, file = "N2_sites.bed", col.names = F, row.names = F, sep = "\t", quote = F)
write.table(N3, file = "N3_sites.bed", col.names = F, row.names = F, sep = "\t", quote = F)