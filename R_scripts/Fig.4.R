libraries <- c("ggplot2", "ggpolypath", "ggpubr", "ggrepel", "data.table", 
               "rstatix", "tximport", "DESeq2", "GenomicFeatures", "vsn", 
               "pheatmap", "RColorBrewer", "EnhancedVolcano", "venn", 
               "ggplotify", "viridis", "plyr")

invisible(lapply(libraries, library, character.only = TRUE))

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

#DESeq2

# Define and set the working directory
wd <- "/Volumes/matthews/Guy/Raw_data/monomorph/data/salmon"
setwd(wd)

# Read in the sample data
samples <- read.table(file.path(wd,"samples.txt"), header=TRUE)

# Set the row names
rownames(samples) <- samples$sample
samples_clonal <- samples[ which(samples$selection=='Clonal') ,]

# Import salmon quants
files <- file.path(wd,"quants", samples_clonal$sample, "quant.sf")
txdb = makeTxDbFromGFF('TriTrypDB-59_TbruceiTREU927.gff')
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
txi <- tximport(files, type="salmon", tx2gene=tx2gene)

names(files) <- samples_clonal$sample
samples_clonal$replicate <- as.factor(samples_clonal$replicate)

ddsTxi_clone <- DESeqDataSetFromTximport(txi,
                                         colData = samples_clonal,
                                         design = ~ stage)

#Keep genes above 10
keep_clone <- rowSums(counts(ddsTxi_clone)) >= 10
dds_clone <- ddsTxi_clone[keep_clone,]

#Run for bulk
dds_clone <- DESeq(dds_clone)
res <- results(dds_clone, name="stage_Pleomorph_vs_Monomorph")
summary(res)

#Run per clone
dds_clone$group <- factor(paste0(dds_clone$clone, dds_clone$stage))
design(dds_clone) <- ~ group
dds_clone <- DESeq(dds_clone)

# Define the clone names
clones <- c("A2", "A4", "A5", "A6", "A7")

# Use lapply to get results for each clone
results_list <- lapply(clones, function(clone) {
  results(dds_clone, contrast=c("group", paste0(clone, "Pleomorph"), paste0(clone, "Monomorph")))
})

# Assign names to the list elements
names(results_list) <- clones

res.df <- setDT(as.data.frame(res), keep.rownames = TRUE)[]
setnames(res.df , old=c("log2FoldChange","padj"), new=c("res_log2FoldChange","res_padj"))
res.df [, c("baseMean","lfcSE", "pvalue", "stat"):=NULL]

df_list <- lapply(clones, function(clone) {
  df <- setDT(as.data.frame(results_list[[clone]]), keep.rownames = TRUE)[]
  setnames(df, old=c("log2FoldChange","padj"), new=c(paste0(clone, "_log2FoldChange"),paste0(clone, "_padj")))
  df[, c("baseMean","lfcSE", "pvalue", "stat"):=NULL]
  df
})

names(df_list) <- clones

combi_df <- Reduce(function(x, y) merge(x, y, by = "rn"), c(list(res.df), df_list))

long_padj <- melt(combi_df, measure.vars = c("A2_padj", "A4_padj", 
                                             "A5_padj", "A6_padj", 
                                             "A7_padj"), id.vars = "rn" )

long_lfc <- melt(combi_df, measure.vars = c("A2_log2FoldChange", "A4_log2FoldChange", 
                                            "A5_log2FoldChange", "A6_log2FoldChange", 
                                            "A7_log2FoldChange"), id.vars = "rn" )

# Define the clone names
clones <- c("res", "A2", "A4", "A5", "A6", "A7")

# Filter significant differences for each clone
sigdif_list <- lapply(clones, function(clone) {
  combi_df[which(combi_df[[paste0(clone, "_padj")]] == 0 | 
                   combi_df[[paste0(clone, "_padj")]] < 0.05 & 
                   combi_df[[paste0(clone, "_log2FoldChange")]] < -1 | 
                   combi_df[[paste0(clone, "_log2FoldChange")]] > 1),]
})

# Assign names to the list elements
names(sigdif_list) <- clones

common_elements <- Reduce(intersect, list(sigdif_list$A2$rn, sigdif_list$A4$rn, sigdif_list$A5$rn, sigdif_list$A6$rn, sigdif_list$A7$rn));
common_elements

#Write combi_df
#write.csv(combi_df, 
#          file="combi_results.csv")


#RBP10 bound enrichment calculation (https://rpubs.com/jrgonzalezISGlobal/enrichment)
#Genes bound by RBP10 and DE in monomorph = 7, Genes bound by RBP10 = 261, Sig. dif in monomorph = 21, Genes analysed = 9551
#
#deTable <-  matrix(c(8, 261, 21, 9551),
#                   nrow = 2,
#                   dimnames = list(DE=c("yes","no"),
#                                   GeneSet=c("in","out")))
#
#fisher.test(deTable, alternative = "greater")

# Transformation and clustering
vsd <- vst(dds_clone, blind=FALSE)
select <- order(rowMeans(counts(dds_clone,normalized=TRUE)),
                decreasing=TRUE)[1:20]

df <- as.data.frame(colData(dds_clone)[,c("stage","clone")])

# Create a distance matrix for pca
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pcaData <- plotPCA(vsd, intgroup=c("stage", "clone"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

b <- ggplot(pcaData, aes(PC1, PC2, color=stage, shape=clone)) +
  geom_point(size=3) +
  ylab(paste0("PC1: ",percentVar[2],"% variance")) +
  xlab(paste0("PC2: ",percentVar[1],"% variance")) + 
  coord_fixed() +
  theme_bw () +
  theme(text = element_text(size = 15)) + 
  annotate("rect", xmin=c(-12,1), xmax=c(-1,9), ymin=c(-5,-5) , ymax=c(7,7), alpha=0, color=c("#F8766D", "#00BFC4"))

# Extract data for common elements heatmap
mat <- assay(vsd)[ common_elements, ]

mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("stage","clone")])

d <- as.ggplot(pheatmap(mat, annotation_col = anno))

# Extract data for QS pathway heatmap
QS <- (c("Tb927.2.1810","Tb927.2.2720","Tb927.2.4020","Tb927.3.4560","Tb927.4.670","Tb927.4.3620","Tb927.4.3630","Tb927.4.3640","Tb927.4.3650","Tb927.5.3580","Tb927.6.2300","Tb927.6.2360","Tb927.7.2100","Tb927.7.7160","Tb927.8.2860","Tb927.9.4080","Tb927.9.7550","Tb927.9.13530","Tb927.10.5930","Tb927.10.5940","Tb927.10.5950","Tb927.10.12100","Tb927.10.15020","Tb927.10.16120","Tb927.11.290","Tb927.11.300","Tb927.11.750","Tb927.11.760","Tb927.11.1640","Tb927.11.2250","Tb927.11.3650","Tb927.11.4610","Tb927.11.6600","Tb927.11.11470","Tb927.11.11480","Tb927.10.12090","Tb927.1.1930","Tb927.11.9270","Tb927.6.4220","Tb927.8.1530","Tb927.10.1740","Tb927.10.1750","Tb927.10.2030","Tb927.10.2530","Tb927.10.12110","Tb927.11.1480","Tb927.11.6610","Tb927.7.2660","Tb927.4.5390","Tb927.8.7020","Tb927.11.2500","Tb927.8.8330","Tb927.11.3570","Tb927.6.400","Tb927.11.6590","Tb927.3.2090","Tb927.3.3410","Tb927.11.12850","Tb927.3.4750","Tb927.10.12260","Tb927.1.2100","Tb927.8.2780"))
mat <- assay(vsd)[ QS, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("stage","isolate")])

sa <- as.ggplot(pheatmap(mat, annotation_col = anno, fontsize_row = 3))

diffs <- list(sigdif_list$A2$rn, sigdif_list$A4$rn, sigdif_list$A5$rn, sigdif_list$A6$rn, sigdif_list$A7$rn)
c <- venn(diffs, snames = c("A2", "A4", "A5", "A6", "A7"), ilabels = "counts", zcolor = c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"), ellipse = T,  ggplot = T, borders = T, box = F, ilcs = 1, sncs = 2)

setwd("/Volumes/matthews/Guy/Raw_data/monomorph/data/salmon")

common_go <- read.csv("common_DE.csv")
common_go_slim <- filter(common_go, Benjamini <= 0.05)

sd <- ggplot(common_go_slim, aes(x=Fold.enrichment, y=reorder(ID, -Fold.enrichment), fill = Benjamini, size = Result.count)) +
  geom_point(shape=21) +
  scale_fill_viridis(discrete=F, option="magma") +
  scale_size(range = c(2, 6), name = "Result count") +
  theme_minimal() + 
  ylab("") +
  xlab("Fold enrichment") +
  theme(text = element_text(size = 18)) 

setwd("/Volumes/matthews/Guy/Raw_data/monomorph/data/figures/Fig.4")
qpcr <- read.csv("rbp10_zc3h20_qpcr.csv")
qpcr_out <- split( qpcr , f = qpcr$clone)

sc <- ggboxplot(qpcr_out$A7, x = "dox", y = "ct",
                color = "dox", palette =c("#00AFBB", "#E7B800"),
                add = "jitter", facet.by = c("gene")) +
  theme(text = element_text(size = 18))

#Clonal selection

#Clones A1-A7

A1_7 <- read.csv("clone_BHI_GC.csv")
A1_7$BHI <- factor(A1_7$BHI, levels=c("HMI-9", "HMI-9:BHI"))
A1_7 <-filter(A1_7, clone != 'A1')
A1_7 <-filter(A1_7, clone != 'A3')
A1_7$Hours <- as.numeric(A1_7$Hours)
A1_7$BHI <- factor(A1_7$BHI, levels=c("HMI-9", "HMI-9:BHI"))
A1_7$Competence <- factor(A1_7$Competence, levels=c("Pleomorph", "Monomorph"))

a <- ggline(A1_7, x = "Hours", y = "Norm_density", add = c("mean_se", "jitter"),
            color = "clone", size = 1.5,
            facet.by = c("Competence", "BHI"), short.panel.labs = T, panel.labs = list(BHI = c("HMI-9", "HMI-9:BHI"))) +
  scale_y_continuous(label=scientific_10, limits = c(0, 3200000)) +
  xlab ("Hours") +
  ylab ("Parasites / ml") +
  theme(text = element_text(size = 15))

rb_zc <- read.csv("RBP10_ZC3H20_BHI_GC.csv")

out <- split( rb_zc , f = rb_zc$Conc )
rb_zc_slim <- rbind(out$'0', out$'0.002')

rb_zc_slim$Conc <- as.character(rb_zc_slim$Conc)
rb_zc_slim$BHI <- factor(rb_zc_slim$BHI, levels=c("HMI-9", "HMI-9:BHI"))

rb_zc_slim %>%
  group_by(Dox, BHI, Hours) %>%
  get_summary_stats(Norm_density, type = "mean_sd")

res.aov <- rb_zc_slim %>% 
  group_by(BHI) %>% 
  anova_test(dv = Norm_density, wid = Flask, within = Hours, between = Dox)

pwc <- rb_zc_slim %>%
  group_by(BHI, Hours, Gene) %>% 
  pairwise_t_test(
    Norm_density ~ Conc, paired = F,
    p.adjust.method = "bonferroni"
  )

pwc <- rb_zc_slim %>%
  group_by(BHI, Hours, Gene) %>%
  pairwise_t_test(Norm_density ~ Dox, paired = F, p.adjust.method = "bonferroni") 

pwc <- pwc %>% add_xy_position(x = "Hours")

e <- ggline(rb_zc_slim, x = "Hours", y = "Norm_density", add = c("mean_se", "jitter"),
            color = "Dox", size = 1, palette = c("black", "#E41A1C"), 
            facet.by = c("Gene", "BHI"), short.panel.labs = T) +
  scale_y_continuous(label=scientific_10, limits = c(0, 3200000)) +
  xlim(c(0,78)) +
  xlab ("Hours") +
  ylab ("Parasites / ml") +
  theme(text = element_text(size = 15)) +
  stat_pvalue_manual(pwc, x = "Hours", label = "p.adj.signif", tip.length = 0, hide.ns = T, size = 5) 

#serial dilution plot
rb_zc$Conc <- as.character(rb_zc$Conc)

rb_zc$Conc <- factor(rb_zc$Conc, levels=c("0.2", "0.02", "0.002", "2e-04", "2e-05", "2e-06", "2e-07", "0"))

sb <- ggline(rb_zc, x = "Hours", y = "Norm_density", add = c("mean_se", "jitter"),
             color = "Conc", size = 1, 
             facet.by = c("Gene", "BHI"), short.panel.labs = T) +
  scale_y_continuous(label=scientific_10, limits = c(0, 3200000)) +
  xlim(c(0,78)) +
  xlab ("Hours") +
  ylab ("Parasites / ml") +
  theme(text = element_text(size = 15)) + 
  scale_color_discrete(name = "ng/ml")

#RBP10 ZC3H20 IFA
rb_zc_ifa <- read.csv("RBP10_ZC3H20_ifa.csv")

rb_zc_ifa$total <- (rb_zc_ifa$X1K1N + rb_zc_ifa$X2K1N + rb_zc_ifa$X1K1N + rb_zc_ifa$other)
rb_zc_ifa$dividing <- ((rb_zc_ifa$X2K2N + rb_zc_ifa$X2K1N) / rb_zc_ifa$total)*100
rb_zc_ifa$pad <- ((rb_zc_ifa$pad) / rb_zc_ifa$total)*100

se <- ggboxplot(rb_zc_ifa, x = 'line', y = 'dividing', facet.by = c("Gene"), short.panel.labs = T) + 
  facet_wrap(c("gene"), scales = 'free_x') +
  theme(text = element_text(size = 15)) +
  ylab("2K1N/2K2N cells (%)") +
  xlab("") +
  ggtitle ("") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))

sf <- ggboxplot(rb_zc_ifa, x = 'line', y = 'pad', facet.by = c("Gene"), short.panel.labs = T) + 
  facet_wrap(c("gene"), scales = 'free_x') +
  theme(text = element_text(size = 15)) +
  ylab("PAD+ (%)") +
  xlab("") +
  ggtitle ("") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))

tiff("/Volumes/matthews/Guy/Raw_data/monomorph/data/figures/Fig.4/Fig.4.tiff", units="in", width=15, height=15, res=300)
ggarrange(a, b, c, d, e, ncol = 2, nrow = 3, common.legend = F, legend="right", align = c("hv"), labels = "auto", font.label = list(size = 14, color = "black", face = "bold", family = NULL))
dev.off()

tiff("/Volumes/matthews/Guy/Raw_data/monomorph/data/figures/Fig.4/Fig.S6.tiff", units="in", width=15, height=15, res=300)
ggarrange(sa, sb, sc, sd, se, sf, ncol = 2, nrow = 3,  common.legend = F, legend="right", align = c("hv"), labels = "auto", font.label = list(size = 14, color = "black", face = "bold", family = NULL))
dev.off()
