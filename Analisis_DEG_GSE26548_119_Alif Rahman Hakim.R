# 1. Panggil Library (Sekarang harusnya sudah bisa)
library(GEOquery)
library(limma)

# 2. Menarik data GSE26548 dari NCBI
# Tunggu sampai bar progres di Console selesai 100%
gse <- getGEO("GSE26548", GSEMatrix = TRUE)

# 3. Mengambil objek utama
gset <- gse[[1]]

# 4. Tampilkan informasi data
gset

# 1. Menyiapkan desain eksperimen (GSE26548 memiliki 3 grup: WT, nsp1, nsp2)
# Masing-masing grup punya 3 ulangan (total 9 sampel)
gset$group <- c("WT", "WT", "WT", "nsp1", "nsp1", "nsp1", "nsp2", "nsp2", "nsp2")
gset$group <- as.factor(gset$group)

# 1. Instal library visualisasi (jika belum ada)
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)

# 2. Membuat Volcano Plot
# Gen dengan logFC > 1 dan P.Value < 0.05 dianggap signifikan
tT$diffexpressed <- "No"
tT$diffexpressed[tT$logFC > 1 & tT$adj.P.Val < 0.05] <- "Up"
tT$diffexpressed[tT$logFC < -1 & tT$adj.P.Val < 0.05] <- "Down"

ggplot(data=tT, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("blue", "black", "red")) +
  labs(title="Volcano Plot: nsp1 vs WT (GSE26548)",
       subtitle="Red: Up-regulated, Blue: Down-regulated")
# 2. Membuat model statistik
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gset$group)

# 3. Analisis Perbandingan (Misal kita bandingkan nsp1 vs WT)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(nsp1_vs_WT = nsp1 - WT, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# 4. Mengambil tabel hasil DEG
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
View(tT) # Untuk melihat tabel gen yang paling signifikan

# Instal pheatmap jika belum ada
if (!require("pheatmap", quietly = TRUE)) install.packages("pheatmap")
library(pheatmap)

# Ambil ekspresi untuk 50 gen teratas
top50_genes <- tT[1:50, ]
exp_matrix <- exprs(gset)[rownames(top50_genes), ]

# Gambar Heatmap
pheatmap(exp_matrix, show_colnames = T, main = "Top 50 DEGs Heatmap",
         annotation_col = data.frame(Group = gset$group, row.names = colnames(exp_matrix)))

# Heatmap yang lebih rapi
pheatmap(exp_matrix, 
         show_colnames = T, 
         main = "Top 50 DEGs Heatmap: nsp1 vs WT",
         cluster_cols = FALSE,        # Menjaga urutan grup tetap rapi
         fontsize_row = 7,            # Mengecilkan tulisan nama gen
         fontsize_col = 8,            # Mengecilkan tulisan nama sampel
         annotation_col = data.frame(Group = gset$group, row.names = colnames(exp_matrix)),
         color = colorRampPalette(c("blue", "white", "red"))(100)) # Warna standar bioinfo

if (!require("enrichR", quietly = TRUE)) install.packages("enrichR")
library(enrichR)

# List database yang tersedia
dbs <- c("GO_Biological_Process_2023", "KEGG_2019_Plants")

# Jalankan enrichment (menggunakan nama gen dari tabel tT kamu)
enriched <- enrichr(rownames(tT[1:100,]), dbs)

# Plot GO
plotEnrich(enriched[[1]], showTerms = 10, numChar = 40, y = "Count", orderBy = "P.value")