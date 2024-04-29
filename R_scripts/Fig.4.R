library(ggplot2)
library(ggpolypath)
library(ggpubr)
library(ggrepel)
library(data.table)
library(rstatix)
library (tximport)
library (DESeq2)
library(GenomicFeatures)
library(vsn)
library(pheatmap)
library(RColorBrewer)
library(data.table)
library (EnhancedVolcano)
library (ggpubr)
library(data.table)
library(venn)
library(ggplotify)
library (viridis)

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

#DESeq2

####################################################################
# Filtered clonal analysis
####################################################################

#Read in the meta data and modify it for tximport
wd <- ("/Users/goldriev/keep/Store/PHD/Monomorph_wet/Clonal_selection/BHI_selection/RNA_seq/salmon")
dir <- system.file("extdata", package="tximportData")
sample_dir <- ("/Users/goldriev/keep/Store/PHD/Monomorph_wet/Clonal_selection/BHI_selection/RNA_seq/salmon")
setwd("/Users/goldriev/keep/Store/PHD/Monomorph_wet/Clonal_selection/BHI_selection/RNA_seq/salmon")
samples <- read.table(file.path(wd,"samples.txt"), header=TRUE)
rownames(samples) <- samples$sample

samples_clonal <- samples[ which(samples$selection=='Clonal') ,]

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
dds_clone$group
#Run per clone
dds_clone$group <- factor(paste0(dds_clone$isolate, dds_clone$stage))
design(dds_clone) <- ~ group
dds_clone <- DESeq(dds_clone)

#Results per clone
A2 <- results(dds_clone, contrast=c("group", "A2Pleomorph", "A2Monomorph"))
A4 <- results(dds_clone, contrast=c("group", "A4Pleomorph", "A4Monomorph"))
A5 <- results(dds_clone, contrast=c("group", "A5Pleomorph", "A5Monomorph"))
A6 <- results(dds_clone, contrast=c("group", "A6Pleomorph", "A6Monomorph"))
A7 <- results(dds_clone, contrast=c("group", "A7Pleomorph", "A7Monomorph"))

res.df <- setDT(as.data.frame(res), keep.rownames = TRUE)[]
setnames(res.df , old=c("log2FoldChange","padj"), new=c("res_log2FoldChange","res_padj"))
res.df [, c("baseMean","lfcSE", "pvalue", "stat"):=NULL]

A2.df <- setDT(as.data.frame(A2), keep.rownames = TRUE)[]
setnames(A2.df , old=c("log2FoldChange","padj"), new=c("A2_log2FoldChange","A2_padj"))
A2.df [, c("baseMean","lfcSE", "pvalue", "stat"):=NULL] 

A4.df <- setDT(as.data.frame(A4), keep.rownames = TRUE)[]
setnames(A4.df , old=c("log2FoldChange","padj"), new=c("A4_log2FoldChange","A4_padj"))
A4.df [, c("baseMean","lfcSE", "pvalue", "stat"):=NULL] 

A5.df <- setDT(as.data.frame(A5), keep.rownames = TRUE)[]
setnames(A5.df , old=c("log2FoldChange","padj"), new=c("A5_log2FoldChange","A5_padj"))
A5.df [, c("baseMean","lfcSE", "pvalue", "stat"):=NULL] 

A6.df <- setDT(as.data.frame(A6), keep.rownames = TRUE)[]
setnames(A6.df , old=c("log2FoldChange","padj"), new=c("A6_log2FoldChange","A6_padj"))
A6.df [, c("baseMean","lfcSE", "pvalue", "stat"):=NULL] 

A7.df <- setDT(as.data.frame(A7), keep.rownames = TRUE)[]
setnames(A7.df , old=c("log2FoldChange","padj"), new=c("A7_log2FoldChange","A7_padj"))
A7.df [, c("baseMean","lfcSE", "pvalue", "stat"):=NULL] 

combi_df <- merge(res.df, A2.df, by = ("rn"))
combi_df <- merge(combi_df, A4.df, by = ("rn"))
combi_df <- merge(combi_df, A5.df, by = ("rn"))
combi_df <- merge(combi_df, A6.df, by = ("rn"))
combi_df <- merge(combi_df, A7.df, by = ("rn"))

long_padj <- melt(combi_df, measure.vars = c("A2_padj", "A4_padj", 
                                             "A5_padj", "A6_padj", 
                                             "A7_padj"), id.vars = "rn" )

long_lfc <- melt(combi_df, measure.vars = c("A2_log2FoldChange", "A4_log2FoldChange", 
                                            "A5_log2FoldChange", "A6_log2FoldChange", 
                                            "A7_log2FoldChange"), id.vars = "rn" )

res_sigdif <- combi_df[ which(combi_df$res_padj == 0 | combi_df$res_padj < 0.05 & combi_df$res_log2FoldChange < -1 | combi_df$res_log2FoldChange > 1),]
A2_sigdif <- combi_df[ which(combi_df$A2_padj == 0 | combi_df$A2_padj < 0.05 & combi_df$A2_log2FoldChange < -1 | combi_df$A2_log2FoldChange > 1),]
A4_sigdif <- combi_df[ which(combi_df$A4_padj == 0 | combi_df$A4_padj < 0.05 & combi_df$A4_log2FoldChange < -1 | combi_df$A4_log2FoldChange > 1),]
A5_sigdif <- combi_df[ which(combi_df$A5_padj == 0 | combi_df$A5_padj < 0.05 & combi_df$A5_log2FoldChange < -1 | combi_df$A5_log2FoldChange > 1),]
A6_sigdif <- combi_df[ which(combi_df$A6_padj == 0 | combi_df$A6_padj < 0.05 & combi_df$A6_log2FoldChange < -1 | combi_df$A6_log2FoldChange > 1),]
A7_sigdif <- combi_df[ which(combi_df$A7_padj == 0 | combi_df$A7_padj < 0.05 & combi_df$A7_log2FoldChange < -1 | combi_df$A7_log2FoldChange > 1),]

#Write combi
#write.csv(combi_df, 
#          file="combi_results.csv")

require(plyr)

common_elements <- Reduce(intersect, list(A2_sigdif$rn, A4_sigdif$rn, A5_sigdif$rn, A6_sigdif$rn, A7_sigdif$rn));
common_elements
#RBP10 bound enrichment calculation (https://rpubs.com/jrgonzalezISGlobal/enrichment)
#Genes bound by RBP10 and DE in monomorph = 7, Genes bound by RBP10 = 261, Sig. dif in monomorph = 21, Genes analysed = 9551

deTable <-  matrix(c(8, 261, 21, 9551),
                   nrow = 2,
                   dimnames = list(DE=c("yes","no"),
                                   GeneSet=c("in","out")))

fisher.test(deTable, alternative = "greater")

# Transformation and clustering

vsd <- vst(dds_clone, blind=FALSE)
rld <- rlog(dds_clone, blind=FALSE)
head(assay(vsd), 3)

ntd <- normTransform(dds_clone)

select <- order(rowMeans(counts(dds_clone,normalized=TRUE)),
                decreasing=TRUE)[1:20]

df <- as.data.frame(colData(dds_clone)[,c("stage","isolate")])

sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pcaData <- plotPCA(vsd, intgroup=c("stage", "isolate"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

b <- ggplot(pcaData, aes(PC1, PC2, color=isolate, shape=stage)) +
  geom_point(size=3) +
  ylab(paste0("PC1: ",percentVar[1],"% variance")) +
  xlab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme(text = element_text(size = 15)) +
  theme_bw ()

mat <- assay(vsd)[ common_elements, ]

mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("stage","isolate")])

d <- as.ggplot(pheatmap(mat, annotation_col = anno))

SIF <- (c("Tb927.2.1810","Tb927.2.2720","Tb927.2.4020","Tb927.3.4560","Tb927.4.670","Tb927.4.3620","Tb927.4.3630","Tb927.4.3640","Tb927.4.3650","Tb927.5.3580","Tb927.6.2300","Tb927.6.2360","Tb927.7.2100","Tb927.7.7160","Tb927.8.2860","Tb927.9.4080","Tb927.9.7550","Tb927.9.13530","Tb927.10.5930","Tb927.10.5940","Tb927.10.5950","Tb927.10.12100","Tb927.10.15020","Tb927.10.16120","Tb927.11.290","Tb927.11.300","Tb927.11.750","Tb927.11.760","Tb927.11.1640","Tb927.11.2250","Tb927.11.3650","Tb927.11.4610","Tb927.11.6600","Tb927.11.11470","Tb927.11.11480","Tb927.10.12090","Tb927.1.1930","Tb927.11.9270","Tb927.6.4220","Tb927.8.1530","Tb927.10.1740","Tb927.10.1750","Tb927.10.2030","Tb927.10.2530","Tb927.10.12110","Tb927.11.1480","Tb927.11.6610","Tb927.7.2660","Tb927.4.5390","Tb927.8.6930","Tb927.9.6100","Tb927.9.6090","Tb927.8.7020","Tb927.11.2500","Tb927.8.8330","Tb927.11.3570","Tb927.6.400","Tb927.11.6590","Tb927.3.2090","Tb927.3.3410","Tb927.11.12850","Tb927.3.4750","Tb927.10.12260","Tb927.1.2100","Tb927.8.2780"))
mat <- assay(vsd)[ SIF, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("stage","isolate")])

sa <- as.ggplot(pheatmap(mat, annotation_col = anno))

diffs <- list(A2_sigdif$rn, A4_sigdif$rn, A5_sigdif$rn, A6_sigdif$rn, A7_sigdif$rn)
c <- venn(diffs, snames = c("A2", "A4", "A5", "A6", "A7"), ilab=T, zcolor = c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"), ellipse = T,  ggplot = T, borders = T, box = F, ilcs = 1, sncs = 2)

setwd("/Volumes/matthews/Guy/Raw_data/monomorph/data/go_term")

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

setwd("/Users/goldriev/Google Drive/My Drive/Developmental_competence_ms/Data/RBP10_ZC3H20/qpcr/")
qpcr <- read.csv("final_rbp10_zc3h20_qpcr.csv")
qpcr_out <- split( qpcr , f = qpcr$clone)

sc <- ggboxplot(qpcr_out$A7, x = "dox", y = "ct",
          color = "dox", palette =c("#00AFBB", "#E7B800"),
          add = "jitter", facet.by = c("gene")) +
  theme(text = element_text(size = 18))


#Clonal selection

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

#Clones A1-A7

setwd("/Users/goldriev/keep/Store/PHD/Monomorph_wet/Clonal_selection/BHI_selection/A1-A7")
A1_7 <- read.csv("BHI_GC.csv")
A1_7$BHI <- factor(A1_7$BHI, levels=c("HMI-9", "HMI-9:BHI"))
A1_7 <-filter(A1_7, Selection != 'A1')
A1_7 <-filter(A1_7, Selection != 'A3')
A1_7$Hours <- as.numeric(A1_7$Hours)
A1_7$BHI <- factor(A1_7$BHI, levels=c("HMI-9", "HMI-9:BHI"))

a <- ggline(A1_7, x = "Hours", y = "Norm_density", add = c("mean_se", "jitter"),
            color = "Selection", size = 1.5,
            facet.by = c("Competence", "BHI"), short.panel.labs = T, panel.labs = list(BHI = c("HMI-9", "HMI-9:BHI"))) +
  scale_y_continuous(label=scientific_10, limits = c(0, 3200000)) +
  xlab ("Hours") +
  ylab ("Parasites / ml") +
  theme(text = element_text(size = 15))

setwd("/Users/goldriev/Google Drive/My Drive/Developmental_competence_ms/Data/RBP10_ZC3H20/")
rb_zc <- read.csv("RBP10_ZC3H20_17.1.24.csv")

out <- split( rb_zc , f = rb_zc$Conc )
rb_zc_slim <- rbind(out$'0', out$'0.002')

rb_zc_slim$Conc <- as.character(rb_zc_slim$Conc)
rb_zc_slim$BHI <- factor(rb_zc_slim$BHI, levels=c("HMI-9", "HMI-9:BHI"))

rb_zc_slim %>%
  group_by(Dox, BHI, Hours) %>%
  get_summary_stats(Norm_density, type = "mean_sd")

rb_zc_out$RBP10 %>%
  group_by(Dox, BHI, Hours) %>%
  identify_outliers(Norm_density)

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

#RBP10 ZC3H20 IFA
rb_zc_ifa <- read.csv("IFA.csv")
rb_zc_ifa <- filter(rb_zc_ifa, gene != "control")

comparisons <- list( c("--", "+-"), c("--", "-+"), c("--", "++"), c("-+", "+-"), c("-+", "++"), c("+-", "++") )

rb_zc_ifa$total <- (rb_zc_ifa$X1K1N + rb_zc_ifa$X2K1N + rb_zc_ifa$X1K1N + rb_zc_ifa$other)
rb_zc_ifa$dividing <- ((rb_zc_ifa$X2K2N + rb_zc_ifa$X2K1N) / rb_zc_ifa$total)*100
rb_zc_ifa$pad <- ((rb_zc_ifa$pad) / rb_zc_ifa$total)*100

f <- ggboxplot(rb_zc_ifa, x = 'line', y = 'dividing', facet.by = c("Gene"), short.panel.labs = T) + 
  facet_wrap(c("gene"), scales = 'free_x') +
  stat_compare_means(comparisons = comparisons, label = "p.signif", method = "t.test",  hide.ns = F) +
  theme(text = element_text(size = 15)) +
  ylab("PAD+ (%)") +
  xlab("") +
  ggtitle ("")


g <- ggboxplot(rb_zc_ifa, x = 'line', y = 'pad', facet.by = c("Gene"), short.panel.labs = T) + 
  facet_wrap(c("gene"), scales = 'free_x') +
  stat_compare_means(comparisons = comparisons, label = "p.signif", method = "t.test",  hide.ns = F) +
  theme(text = element_text(size = 15)) +
  ylab("PAD+ (%)") +
  xlab("") +
  ggtitle ("")

#serial dilution plot
setwd("/Users/goldriev/Google Drive/My Drive/Developmental_competence_ms/Data/RBP10_ZC3H20/")
rb_zc <- read.csv("RBP10_ZC3H20_17.1.24.csv")
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

tiff("/Users/goldriev/Google Drive/My Drive/Developmental_competence_ms/draft_ms/figures/Fig.4/Fig.4.tiff", units="in", width=15, height=15, res=300)
ggarrange(a, b, c, d, e, f, g, ncol = 2, nrow = 4, common.legend = F, legend="right", align = c("hv"), labels = "auto", font.label = list(size = 14, color = "black", face = "bold", family = NULL))
dev.off()

tiff("/Users/goldriev/Google Drive/My Drive/Developmental_competence_ms/draft_ms/figures/Fig.4/Fig.S5.tiff", units="in", width=15, height=15, res=300)
ggarrange(sa, sb, sc, sd, ncol = 2, nrow = 2,  common.legend = F, legend="right", align = c("hv"), labels = "auto", font.label = list(size = 14, color = "black", face = "bold", family = NULL))
dev.off()
