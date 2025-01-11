  Sys.setenv(http_proxy = "http://proxy_address:port", https_proxy = "http://proxy_address:port")
install.packages("BiocManager")
BiocManager::install("GEOquery", force = TRUE)
library(GEOquery)
options(download.file.method = "wininet")  # Windows için
# veya
options(download.file.method = "wget")
Sys.getenv("PATH")
Sys.setenv(PATH = paste("C:\\Rtools\\bin", Sys.getenv("PATH"), sep = ";"))
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("GEOquery")
capabilities("http/ftp")
Sys.unsetenv("http_proxy")
Sys.unsetenv("https_proxy")
options(download.file.method = "wininet")
install.packages("BiocManager")
library(BiocManager)
# GEOquery paketini yükleme
BiocManager::install("GEOquery")
Sys.setenv(http_proxy = "http://proxy.example.com:8080")
Sys.setenv(https_proxy = "http://proxy.example.com:8080")
# Gerekli paketleri yükleyin
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("GEOquery")
install.packages("BiocManager")
BiocManager::install(version = "3.16") # Geçerli Bioconductor sürümünüzü kontrol edin
BiocManager::install(version = "3.16") # Geçerli Bioconductor sürümünüzü kontrol edin
library(caret)
library(pheatmap)
# Isı haritası
pheatmap(exprs_data[1:50, ], cluster_rows = TRUE, cluster_cols = TRUE)
library(clusterProfiler)
library(pROC)
library(limma)
# Tasarım matrisi
group <- factor(c("control", "treated", "control", "treated"))  # Örnek bilgilerini düzenleyin
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
# Farklı ifade edilen genleri analiz edin
fit <- lmFit(exprs_data, design)
library(limma)
# Tasarım matrisi
group <- factor(c("control", "treated", "control", "treated"))  # Örnek bilgilerini düzenleyin
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
# Farklı ifade edilen genleri analiz edin
fit <- lmFit(exprs_data, design)
rm(list=ls()) #clear extra stuff
setwd("~/Dropbox/Choughs/Plots/")
rm(list=ls()) #clear extra stuff
setwd("C:\\Users\\ESHA\\Documents\\Thesis\\For R")
# Example data
urban_abundance <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)  # Mean abundance of predators in urban habitats
natural_abundance <- c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)  # Mean abundance of predators in natural habitats
# Perform paired t-test
t_test_result <- t.test(urban_abundance, natural_abundance)
# Print the result
print(t_test_result)
rm(list=ls()) #clear extra stuff
setwd("C:\\Users\\ESHA\\Documents\\Thesis\\For R")
library(coda)
library(rjags)
# Auto detect text files and perform LF normalization
* text=auto
# Auto detect text files and perform LF normalization
text=auto
###########################################
#' Check for packages and if necessary install into library
#+ message = FALSE
pkgs <- c("randomForest","data.table", "dplyr", "lubridate","lme4","here")
# Check if packages are already installed
inst <- pkgs %in% installed.packages()
# Install missing packages
if (any(!inst)) install.packages(pkgs[!inst])
# Load packages
pkg_out <- lapply(pkgs, require, character.only = TRUE)
# Set the working directory for the "here" package
here::i_am("EPHI_paper.Rproj")
###########################################
###########################################
#' Check for packages and if necessary install into library
#+ message = FALSE
pkgs <- c("randomForest","data.table", "dplyr", "lubridate","lme4","R2jags","here")
# Check if packages are already installed
inst <- pkgs %in% installed.packages()
# Install missing packages
if (any(!inst)) install.packages(pkgs[!inst])
# Load packages
pkg_out <- lapply(pkgs, require, character.only = TRUE)
# Set the working directory for the "here" package
here::i_am("EPHI_paper.Rproj")
###########################################
library(data.table)
library(dplyr)
library(R2jags)
update.packages(ask = FALSE)
library(R2jags)
pkgs <- c("randomForest", "data.table", "dplyr", "lubridate", "lme4", "R2jags", "here")
inst <- pkgs %in% installed.packages()
if (any(!inst)) install.packages(pkgs[!inst])
pkgs <- c("randomForest", "data.table", "dplyr", "lubridate", "lme4", "R2jags", "here")
inst <- pkgs %in% installed.packages()
if (any(!inst)) install.packages(pkgs[!inst])
lapply(pkgs, library, character.only = TRUE)
library(R2jags)  # JAGS yazılımı kuruluysa hata almazsınız
C:\Program Files\JAGS\JAGS-4.x.y
library(R2jags)
remove.packages("rjags")
remove.packages("R2jags")
install.packages("rjags")
install.packages("R2jags")
###########################################
#' Check for packages and if necessary install into library
#+ message = FALSE
pkgs <- c("randomForest","data.table", "dplyr", "lubridate","lme4","here")
# Check if packages are already installed
inst <- pkgs %in% installed.packages()
# Install missing packages
if (any(!inst)) install.packages(pkgs[!inst])
# Load packages
pkg_out <- lapply(pkgs, require, character.only = TRUE)
# Set the working directory for the "here" package
here::i_am("EPHI_paper.Rproj")
library(here)
here::i_am("EPHI_paper.Rproj")
setwd("C:/Users/User/Documents")
###########################################
###########################################
#' Check for packages and if necessary install into library
#+ message = FALSE
pkgs <- c("randomForest","data.table", "dplyr", "lubridate","lme4","R2jags","here")
# Check if packages are already installed
inst <- pkgs %in% installed.packages()
# Install missing packages
if (any(!inst)) install.packages(pkgs[!inst])
# Load packages
pkg_out <- lapply(pkgs, require, character.only = TRUE)
# Set the working directory for the "here" package
here::i_am("EPHI_paper.Rproj")
# Kurulu olan tüm paketleri listeleme
installed_packages <- installed.packages()
installed_packages_df <- as.data.frame(installed_packages)
# Paket adlarını görmek için
print(installed_packages_df$Package)
remove.packages("rjags")
remove.packages("R2jags")
install.packages("rjags")
install.packages("R2jags")
# Kurulu olan tüm paketleri listeleme
installed_packages <- installed.packages()
installed_packages_df <- as.data.frame(installed_packages)
# Paket adlarını görmek için
print(installed_packages_df$Package)
###########################################
###########################################
#' Check for packages and install them if necessary
#+ message = FALSE
pkgs <- c("randomForest","data.table", "dplyr", "lubridate","lme4","R2jags",
"mcmcOutput","mcmcplots","MCMCvis","here")
# Check if packages are already installed
inst <- pkgs %in% installed.packages()
# Install missing packages
if (any(!inst)) install.packages(pkgs[!inst])
# Load packages
pkg_out <- lapply(pkgs, require, character.only = TRUE)
# Set the working directory for the "here" package
here::i_am("EPHI_paper.Rproj")
###########################################
###########################################
#' Check for packages and if necessary install into library
#+ message = FALSE
pkgs <- c("randomForest","data.table", "dplyr", "lubridate","lme4","R2jags","mcmcOutput","mcmcplots","MCMCvis",
"pastecs","ggplot2","cowplot","gridExtra","scales","here")
# Check if packages are already installed
inst <- pkgs %in% installed.packages()
# Install missing packages
if (any(!inst)) install.packages(pkgs[!inst])
# Load packages
pkg_out <- lapply(pkgs, require, character.only = TRUE)
# Set the working directory for the "here" package
here::i_am("EPHI_paper.Rproj")
integrated_seurat <- FindClusters(integrated_seurat, resolution = 0.6, cluster.name = "rpca_clusters")
"Seurat" %in% installed.packages()
install.packages("Seurat")
library(Seurat)
integrated_seurat <- FindClusters(integrated_seurat, resolution = 0.6, cluster.name = "rpca_clusters")
packageVersion("Seurat")
install.packages("Seurat")
Seurat::FindClusters(integrated_seurat, resolution = 0.6, cluster.name = "rpca_clusters")
class(integrated_seurat)
integrated_seurat <- CreateSeuratObject(counts = raw_data)
install.packages("Seurat")
packageVersion("Seurat")
# 1. Veri yükleme (örnek veri seti)
library(Seurat)
pbmc_data <- Read10X(data.dir = "path_to_your_data_directory")
class(integrated_seurat)
/data_directory/
pbmc_data <- Read10X(data.dir = "C:/Users/User/Documents/data_directory")
integrated_seurat <- FindClusters(integrated_seurat, resolution = 0.6, cluster.name = "rpca_clusters")
>1 dna:chromosome chromosome:Galgal4:1:48849500:49020500:1
1 dna:chromosome chromosome:Galgal4:1:48849500:49020500:1
TGAAGCTTATCGTGTACTCTGTCTTCAGAAGCACAAGGACAATTAACATTCAATAGCAAG
import random
import random
# Import the random module
import random
def fitness(x):
def generate_population(size, x_min, x_max):
# BiocManager'ı yükleyin (eğer kurulu değilse)
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
# Bioconductor ile bir paket yüklemek için:
BiocManager::install("genefilter")  # Örnek paket
a
library(GEOquery)
# Örnek GSE (GEO Series) verisini indirme
gse <- getGEO("GSE12592", GSEMatrix = TRUE)
expression_set <- gse[[1]]  # İlk ExpressionSet nesnesini seçin
# Gen ifadelerini ve meta bilgileri kontrol etme
exprs_data <- exprs(expression_set)  # Gen ifadeleri matrisi
sample_info <- pData(expression_set)  # Örnek bilgileri
feature_info <- fData(expression_set)  # Genetik özellik bilgileri
# İlk birkaç satırı görüntüleme
head(exprs_data)
head(sample_info)
BiocManager::install("limma")  # Eğer limma kurulu değilse
BiocManager::install("clusterProfiler")
library(ggplot2)
# PCA Analizi
pca <- prcomp(t(exprs_data), scale. = TRUE)
# PCA Görselleştirme
pca_data <- data.frame(pca$x)
pca_data$group <- sample_info$group
ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) +
geom_point(size = 4) +
labs(title = "PCA Analysis", x = "PC1", y = "PC2") +
theme_minimal()
BiocManager::install("pheatmap")
BiocManager::install("DESeq2")
# 1. Paketlerin yükleneceği dizini ayarla
.libPaths("C:/Users/User/R/library")
# 2. Bioconductor paketlerini yükle
if (!requireNamespace("BiocManager", quietly = TRUE)) {
install.packages("BiocManager", lib = "C:/Users/User/R/library")
}
BiocManager::install(c("limma", "DESeq2", "clusterProfiler", "pheatmap"), lib = "C:/Users/User/R/library")
# 3. Kütüphaneleri yükle
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
# PCA analizi
pca <- prcomp(t(exprs_data), scale. = TRUE)
pca_data <- data.frame(pca$x)
pca_data$group <- sample_info$group
# PCA görselleştirme
ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) +
geom_point(size = 4) +
labs(title = "PCA Analizi", x = "PC1", y = "PC2") +
theme_minimal()
library(ggplot2)
# Önemli genler
deg_genes$Significant <- ifelse(deg_genes$adj.P.Val < 0.05 & abs(deg_genes$logFC) > 1, "Significant", "Not Significant")
library(clusterProfiler)
# KEGG analizi
kegg_enrichment <- enrichKEGG(gene = gene_list, organism = "hsa", keyType = "kegg")
library(enrichplot)
colnames(sample_info)
sample_info$group <- sample_info$your_group_column  # Doğru sütunu belirtin
sample_info$group <- sample_info$geo_accession  # Doğru sütunu belirtin
design <- model.matrix(~ group, data = sample_info)
fit <- lmFit(exprs_data, design)
pca_data$group <- sample_info$group
# Paket Yükleme ve Yolları Belirleme
.libPaths("C:/Users/User/R/library")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
install.packages("BiocManager")
}
BiocManager::install(c("limma", "pheatmap", "clusterProfiler", "GEOquery"))
# Kütüphaneleri Yükleme
library(GEOquery)
library(limma)
table(sample_info$source_name_ch1)
table(sample_info$characteristics_ch1)
table(sample_info$title)
sample_info$group <- sample_info$source_name_ch1
library(DESeq2)
installed_packages_df((DESeq2))
install.packages(DESeq2)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
install.packages("BiocManager")
}
BiocManager::install("DESeq2")
library(DESeq2)
BiocManager::install("DESeq2", lib = "C:/Users/User/R/library")
R.version.string
installed.packages()
installed.packages(DESeq2)
"DESeq2" %in% rownames(installed.packages())
if (!requireNamespace("BiocManager", quietly = TRUE)) {
install.packages("BiocManager")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
install.packages("BiocManager")}
if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
# Bioconductor kurulumunu kontrol et
if (!requireNamespace("BiocManager", quietly = TRUE)) {
install.packages("BiocManager")
}
# DESeq2 yükle
BiocManager::install("DESeq2", dependencies = TRUE)
# Paket yükleme kontrolü
library(DESeq2)
library(DESeq2)
# Simulated expression data and sample metadata
count_data <- matrix(rpois(1000, lambda = 10), nrow = 100, ncol = 10)
col_data <- data.frame(
condition = factor(rep(c("control", "treatment"), each = 5))
)
rownames(count_data) <- paste0("gene", 1:100)
rownames(col_data) <- colnames(count_data) <- paste0("sample", 1:10)
# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ condition)
# Perform differential expression analysis
dds <- DESeq(dds)
results <- results(dds)
# Filter significant genes
significant_genes <- results[results$padj < 0.05 & abs(results$log2FoldChange) > 1, ]
head(significant_genes)
library(ggplot2)
# Perform PCA
pca <- prcomp(t(count_data), scale. = TRUE)
pca_data <- data.frame(pca$x, condition = col_data$condition)
# Plot PCA
ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
geom_point(size = 4) +
labs(title = "PCA of Gene Expression", x = "PC1", y = "PC2") +
theme_minimal()
library(clusterProfiler)
library(org.Hs.eg.db)
# Example gene list
gene_list <- rownames(significant_genes)
# GO enrichment analysis
go_results <- enrichGO(
gene = gene_list,
OrgDb = org.Hs.eg.db,
keyType = "SYMBOL",
ont = "BP", # Biological Process
pAdjustMethod = "fdr",
qvalueCutoff = 0.05
)
# Visualize results
dotplot(go_results)
library(qqman)
# Giriş gen listesini inceleyin
print(gene_list)
# Örnek genler: CECR2, PIWIL1, CFAP73
print(head(deg_genes))
# Daha az katı kriterler
deg_genes <- results[results$padj < 0.1 & abs(results$log2FoldChange) > 0.5, ]
# DEG kontrolü
if (nrow(deg_genes) == 0) {
print("No significant genes found. Consider relaxing the thresholds.")
} else {
gene_list <- rownames(deg_genes)
}
# Örnek gen listesi
gene_list <- c("CECR2", "PIWIL1", "CFAP73", "SLC9A8", "TXNRD3", "CATSPERD")
# GO zenginleştirme analizi
go_results <- enrichGO(
gene = gene_list,
OrgDb = org.Hs.eg.db,
keyType = "SYMBOL",
ont = "BP",
pAdjustMethod = "fdr",
qvalueCutoff = 0.05
)
# Görselleştirme
if (!is.null(go_results)) {
dotplot(go_results, showCategory = 10) +
ggtitle("GO Enrichment Analysis")
} else {
print("No enrichment results.")
}
# DEG sonuçlarından manuel gen seçimi (örnek)
deg_genes <- results[1:10, ]  # İlk 10 gen
gene_list <- rownames(deg_genes)
# enrichGO tekrar çalıştırma
go_results <- enrichGO(
gene = gene_list,
OrgDb = org.Hs.eg.db,
keyType = "SYMBOL",
ont = "BP",
pAdjustMethod = "fdr",
qvalueCutoff = 0.05
)
deg_genes <- results[results$padj < 0.2 & abs(results$log2FoldChange) > 0.2, ]
if (nrow(deg_genes) == 0) {
print("No significant genes found even with relaxed thresholds.")
} else {
gene_list <- rownames(deg_genes)
print(head(deg_genes))
}
# Örnek bilgilerini kontrol edin
print(head(sample_info))
# İfade matrisini kontrol edin
print(dim(exprs_data))  # Satır: genler, Sütun: örnekler
print(head(exprs_data))
# Grup bilgilerini kontrol edin
print(unique(sample_info$group))
table(sample_info$group)
sample_info$group <- ifelse(grepl("control", sample_info$source_name_ch1), "control", "treatment")
# Örnek gen listesi
gene_list <- c("CECR2", "PIWIL1", "CFAP73", "SLC9A8", "TXNRD3", "CATSPERD")
# enrichGO Analizi
go_results <- enrichGO(
gene = gene_list,
OrgDb = org.Hs.eg.db,
keyType = "SYMBOL",
ont = "BP", # Biological Process
pAdjustMethod = "fdr",
qvalueCutoff = 0.05
)
# Sonuçları Görselleştirme
if (!is.null(go_results)) {
dotplot(go_results, showCategory = 10) +
ggtitle("GO Enrichment Analysis")
} else {
print("No enrichment results.")
}
deg_genes <- results[results$padj < 0.2 & abs(results$log2FoldChange) > 0.2, ]
if (nrow(deg_genes) == 0) {
print("No significant genes found even with relaxed thresholds.")
} else {
gene_list <- rownames(deg_genes)
print(head(deg_genes))
}
set.seed(123)
# Simulating genotype matrix (100 individuals, 10 SNPs)
genotype_matrix <- matrix(sample(0:2, 1000, replace = TRUE), nrow = 100, ncol = 10)
# Simulate phenotype (binary trait)
phenotype <- rbinom(100, 1, prob = 0.5)
# Association test using logistic regression
association_results <- apply(genotype_matrix, 2, function(snp) {
summary(glm(phenotype ~ snp, family = "binomial"))$coefficients[2, 4]  # P-value
})
# Print association results
association_results
library(igraph)
# Compute gene-gene correlation matrix
cor_matrix <- cor(t(count_data), method = "pearson")
# Threshold correlations to create a network
adj_matrix <- (abs(cor_matrix) > 0.7) * 1  # Binary adjacency matrix
# Create igraph object
network <- graph.adjacency(adj_matrix, mode = "undirected", diag = FALSE)
# Plot the network
plot(network, vertex.size = 5, vertex.label = NA, edge.color = "gray")
library(pheatmap)
# Scale the data for better visualization
scaled_data <- t(scale(t(count_data)))
# Create a heatmap
pheatmap(scaled_data,
cluster_rows = TRUE,
cluster_cols = TRUE,
show_rownames = TRUE,
show_colnames = TRUE,
main = "Heatmap of Gene Expression")
library(qqman)
# Veri boyutlarını kontrol edin
print(dim(count_data))  # Gen sayısı ve örnek sayısı
# Örnek veri ile çalışıyorsanız, skalama işleminin uygun olup olmadığını kontrol edin
print(head(scaled_data))
# Eğer herhangi bir hata varsa, count_data matrisinde eksik veya hatalı değerler olabilir.
if (!requireNamespace("qqman", quietly = TRUE)) {
install.packages("qqman")
}
library(qqman)
install.packages("qqman", lib = "C:/Users/User/AppData/Local/R/win-library/4.3")
# Simulated GWAS results
gwas_data <- data.frame(
SNP = paste0("rs", 1:1000),
CHR = rep(1:22, length.out = 1000),
BP = seq(1, 1e6, length.out = 1000),
P = runif(1000, 1e-6, 1)
)
# Manhattan plot
manhattan(gwas_data, main = "GWAS Manhattan Plot", col = c("blue4", "orange3"))
# QQ plot
qq(gwas_data$P, main = "GWAS QQ Plot")
deg_genes <- results[results$padj < 0.2 & abs(results$log2FoldChange) > 0.3, ]
library(clusterProfiler)
library(org.Hs.eg.db)
# Example gene list
gene_list <- rownames(significant_genes)
# GO enrichment analysis
go_results <- enrichGO(
gene = gene_list,
OrgDb = org.Hs.eg.db,
keyType = "SYMBOL",
ont = "BP", # Biological Process
pAdjustMethod = "fdr",
qvalueCutoff = 0.05
)
# Visualize results
dotplot(go_results)
q()
load("~/.RData")
load("~/.RData")
load("~/.RData")
