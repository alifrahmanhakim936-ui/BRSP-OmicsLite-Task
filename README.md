# Analisis Transkriptomik Respon Kekeringan Medicago truncatula

Repository ini berisi analisis fungsional gen pada bintil akar *Medicago truncatula* yang mengalami cekaman kekeringan selama 4 hari (T4) dibandingkan dengan kontrol (T0).

## Deskripsi Proyek
Proyek ini mengidentifikasi Differentially Expressed Genes (DEGs) menggunakan data transkriptomik dari GSE126986. Analisis dilakukan menggunakan bahasa pemrograman R untuk memahami mekanisme pertahanan tanaman legum terhadap stres abiotik.

## Hasil Utama
- **Total DEGs:** 2488 gen (FDR < 0.05, |log2FC| > 1)
- **Up-regulated:** 884 gen
- **Down-regulated:** 1604 gen

### Visualisasi
![Volcano Plot](results/Volcano_Plot_Medicago.png)

## Struktur File
- `scripts/analysis.R`: Script utama pemrosesan data.
- `data/`: Dataset mentah (CSV).
- `results/`: Output visualisasi (Volcano Plot & Tabel).

## Referensi
Sańko-Sawczenko, I., et al. (2019). International Journal of Molecular Sciences.
