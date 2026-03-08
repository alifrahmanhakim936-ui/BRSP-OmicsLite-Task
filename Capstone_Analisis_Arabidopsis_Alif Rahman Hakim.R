# 1. Baca file
data_mt <- read.csv("Supplementary Data Set 1.csv")

# 2. Filter gen yang signifikan (FDR < 0.05 dan logFC > 1 atau < -1)
# Kita bandingkan T4 (Kekeringan) vs T0 (Kontrol)
deg_signifikan <- data_mt[data_mt$T4vsT0...FDR < 0.05 & abs(data_mt$T4vsT0...logFC) > 1, ]

# 3. Print hasilnya (Ini yang nanti kamu tunjukkan di video!)
print(paste("Jumlah gen yang berubah signifikan adalah:", nrow(deg_signifikan)))

# 4. Buat tabel 10 gen teratas untuk bahan presentasi
top_genes <- head(deg_signifikan[order(deg_signifikan$T4vsT0...FDR), ], 10)
print(top_genes[, c("Mt_GeneID", "Mt_Description", "T4vsT0...logFC")])

# 5. Memfilter gen dengan FDR < 0.05 dan |logFC| > 1
# Perbaikan: Menambahkan koma di akhir dan perbandingan > 1
deg_signifikan <- data_mt[data_mt$T4vsT0...FDR < 0.05 & abs(data_mt$T4vsT0...logFC) > 1, ]

# Melihat jumlah gen
print(paste("Jumlah gen signifikan:", nrow(deg_signifikan)))

# 6. Mengambil 10 Gen paling signifikan
# Perbaikan: Nama variabel deg_signifikan sudah benar (pakai 'i')
top_10 <- head(deg_signifikan[order(deg_signifikan$T4vsT0...FDR), ], 10)
print(top_10[, c("Mt_GeneID", "Mt_Description", "T4vsT0...logFC")])

# 7. Load library ggplot2
library(ggplot2)

# Menambah kolom kategori untuk warna di grafik
data_mt$category <- "Not Significant"
data_mt$category[data_mt$T4vsT0...logFC > 1 & data_mt$T4vsT0...FDR < 0.05] <- "Up-regulated"
data_mt$category[data_mt$T4vsT0...logFC < -1 & data_mt$T4vsT0...FDR < 0.05] <- "Down-regulated"

# Membuat Volcano Plot (Perbaikan: Kurung di labs dan typo theme)
ggplot(data_mt, aes(x = T4vsT0...logFC, y = -log10(T4vsT0...PValue), color = category)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c("Down-regulated" = "blue", "Not Significant" = "grey", "Up-regulated" = "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Respon Medicago truncatula (Kekeringan T4 vs T0)",
       x = "Log2 Fold Change",
       y = "-Log10 P-Value") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")

# 1. Pastikan saya ambil 10 gen teratas berdasarkan FDR terkecil
top_10_table <- head(deg_signifikan[order(deg_signifikan$T4vsT0...FDR), ], 10)

# 2. Pilih kolom yang penting saja agar muat di layar
# Kita ambil ID, Deskripsi singkat, dan nilai logFC-nya
clean_table <- top_10_table[, c("Mt_GeneID", "Mt_Description", "T4vsT0...logFC")]

# 3. Ubah nama kolom agar lebih rapi saat dipresentasikan
colnames(clean_table) <- c("Gene ID", "Description", "Log2 Fold Change")

# 4. Munculkan di panel Plot (berjajaran dengan Volcano Plot)
grid.table(clean_table)