library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library (ggpubr)
library(dplyr)
library(ggridges)
library(splitstackshape)
library(tidyr)
library(circlize)
library(ComplexHeatmap)
library(viridis)
library(ggplotify)
library(gridBase)
library(ape)
library(Biostrings)
library(ggtree)
library(treeio)
library(phytools)
library(cowplot)
library(gridExtra)

setwd("/Volumes/matthews/Guy/Raw_data/monomorph/data/circos")

chr <- read.table('chr.txt')
eva_variants <- read.table('T.b.evansitypeA.vcf.variants.shuf', header=F)
evb_variants <- read.table('T.b.evansitypeB.vcf.variants.shuf', header=F)
evi_variants <- read.table('T.b.evansitypeC.vcf.variants.shuf', header=F)
ovi_variants <- read.table('T.b.equitypeOVI.vcf.variants.shuf', header=F)
botat_variants <- read.table('T.b.equitypeBOTAT.vcf.variants.shuf', header=F)

eva_roh <- read.table('evansi_a_roh.txt', header=F)
evb_roh <- read.table('evansi_b_roh.txt', header=F)
evi_roh <- read.table('evansi_c_roh.txt', header=F)
ovi_roh <- read.table('equi_ovi_roh.txt', header=T, col.names = c('chr','start','stop'))
botat_roh <- read.table('equi_botat_roh.txt', header=F)

eva_cnv <- read.table('evansi_a_cnv.txt', header=F)
evb_cnv <- read.table('evansi_b_cnv.txt', header=F)
evi_cnv <- read.table('evansi_c_cnv.txt', header=F)
ovi_cnv <- read.table('equi_ovi_cnv.txt', header=F)
botat_cnv <- read.table('equi_botat_cnv.txt', header=F)

eva_dnds <- read.table('evansi_a_dnds.txt', header=F)
evb_dnds <- read.table('evansi_b_dnds.txt', header=F)
evi_dnds <- read.table('evansi_c_dnds.txt', header=F)
ovi_dnds <- read.table('equi_ovi_dnds.txt', header=F)
botat_dnds <- read.table('equi_botat_dnds.txt', header=F)

variant_list <- list(eva_variants, evb_variants, evi_variants, ovi_variants, botat_variants)
roh_list <- list(eva_roh, evb_roh, evi_roh, ovi_roh, botat_roh)
cnv_list <- list(eva_cnv, evb_cnv, evi_cnv, ovi_cnv, botat_cnv)
dnds_list <- list(eva_dnds, evb_dnds, evi_dnds, ovi_dnds, botat_dnds)

circlize_plot = function() {
  circos.par("track.height" = 0.15, cell.padding = c(0, 0, 0, 0), "clock.wise" = T, "gap.after" = 10)
  circos.par(gap.after = c(rep(1, 10), 20))
  circos.genomicInitialize(chr,sector.names = c("1","2", "3", "4", "5", "6", "7", "8", "9", "10", "11"), labels.cex = 0.8)
  col_mat <-  c("darkred", "darkgreen", "#643F95", "gold", "darkorange")
  circos.genomicDensity(variant_list, col = col_mat, baseline = 0, area = F, window.size = 200000, count_by = c("percent"))
  f = colorRamp2(breaks = c(1.9, 2.1), colors = c("darkblue", "red"))
  
  for (x in cnv_list) {
    circos.genomicTrack(x, ylim =c(0,0.01), track.height = 0.02, panel.fun = function(region, value, ...) {
      i = getI(...)
      circos.genomicRect(region, value, col = f(value[[1]]), border = NA)})}
  
  circos.genomicDensity(dnds_list, col = col_mat, baseline = 0, area = F, window.size = 200000, count_by = c("percent"))
  
  for (x in roh_list) {
    circos.genomicTrack(x, ylim =c(0,0.01), track.height = 0.02, panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = "black", border = NA)
    })}
  draw.sector(0, 360, rou1 = 0.015, center = c(0.705, 0.02), border = NA, lwd = 2, lty = 1, col = "darkred")
  draw.sector(0, 360, rou1 = 0.015, center = c(0.665, 0.02), border = NA, lwd = 2, lty = 1, col = "darkgreen")
  draw.sector(0, 360, rou1 = 0.015, center = c(0.625, 0.02), border = NA, lwd = 2, lty = 1, col = "#643F95")
  draw.sector(0, 360, rou1 = 0.015, center = c(0.585, 0.02), border = NA, lwd = 2, lty = 1, col = "gold")
  draw.sector(0, 360, rou1 = 0.015, center = c(0.545, 0.02), border = NA, lwd = 2, lty = 1, col = "darkorange")
  draw.sector(0, 360, rou1 = 0.015, center = c(0.335, 0.02), border = NA, lwd = 2, lty = 1, col = "darkred")
  draw.sector(0, 360, rou1 = 0.015, center = c(0.295, 0.02), border = NA, lwd = 2, lty = 1, col = "darkgreen")
  draw.sector(0, 360, rou1 = 0.015, center = c(0.255, 0.02), border = NA, lwd = 2, lty = 1, col = "#643F95")
  draw.sector(0, 360, rou1 = 0.015, center = c(0.215, 0.02), border = NA, lwd = 2, lty = 1, col = "gold")
  draw.sector(0, 360, rou1 = 0.015, center = c(0.175, 0.02), border = NA, lwd = 2, lty = 1, col = "darkorange")
  circos.clear()
}

circlize_plot()

text(.27, .06, substitute(paste(bold("ROH"))), family = "helvetica", cex = 0.6)
text(.43, .06, substitute(paste(bold("dN/dS"))), family = "helvetica", cex = 0.6)
text(.63, .06, substitute(paste(bold("CNV"))), family = "helvetica", cex = 0.6)
text(.81, .06, substitute(paste(bold("SNP"))), family = "helvetica", cex = 0.6)

p1 <- recordPlot() 
ggdraw(p1)

tree <- treeio::read.newick("~/Google Drive/My Drive/Developmental_competence_ms/Data/RpPVmHuFdPEBY8WfDi-Gww_newick.txt", node.label="support")
root <- rootnode(tree)  
#tree <- treeio::read.newick("/Volumes/matthews/Guy/Raw_data/monomorph/data/tree/strict_snps.genotyped.fasta.varsites.phy.contree", node.label="support")

a <- ggtree(tree, layout="circular") + 
  geom_tiplab(size = 1.5) +
  geom_point2(aes(subset=!isTip & node != root & as.numeric(support) < 100, 
                  fill=cut(support, c(0, 50, 75, 100))), 
              shape=21, size=2) + 
  geom_hilight(node=c(2,3), fill='orange', alpha=0.3, extend=0.1) +
  geom_hilight(node=c(5,6,7,8), fill='gold', alpha=0.3, extend=0.1) +
  geom_hilight(node=c(23,24), fill="darkred", alpha=0.5, extend=0.085) +
  geom_hilight(node=c(22), fill="darkred", alpha=0.5, extend = 0.07) +
  geom_hilight(node=c(56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83), 
               fill="darkgreen", alpha=0.3, extend=0.075) +
  theme_tree(legend.position=c(0.5, 0.4)) + 
  scale_fill_manual(values=c("white", "grey", "black"), guide='legend', 
                    name='Bootstrap Percentage', 
                    breaks=c('(75,100]', '(50,75]', '(25,50]'), 
                    labels=expression("90-100", "70-90", "0-0")) +
  guides(fill = F) + 
  geom_treescale(offset = -2, width = 0.05, linesize = 1, x = 0, y = 0)

fig1_plots = list(p1,a)

tiff("/Users/goldriev/Google Drive/My Drive/Developmental_competence_ms/draft_ms/figures/Fig.1/Fig.1.tiff", units="in", width=6, height=10, res=300)
ggarrange(plotlist = fig1_plots, ncol = 1, nrow = 2, labels = c("b", "a"), font.label = list(size = 14, color = "black", face = "bold", family = NULL))
dev.off()

setwd("/Volumes/matthews/Guy/Raw_data/monomorph/data/go_term/")
ovi_cnv_go <- read.csv("ovi_cnv_go.csv")
botat_roh_go <- read.csv("botat_roh_go.csv")
ovi_roh_go <- read.csv("ovi_roh_go.csv")
common_go <- read.csv("common_dnds_go.csv")
ev_a_go <- read.csv("evansi_a_dnds_go.csv")
ev_b_go <- read.csv("evansi_b_dnds_go.csv")
ev_c_go <- read.csv("evansi_c_dnds_go.csv")
ovi_go <- read.csv("equi_ovi_dnds_go.csv")
botat_go <- read.csv("equi_botat_dnds_go.csv")

ovi_cnv_go_slim <- filter(ovi_cnv_go, Benjamini <= 0.05)
ovi_roh_go_slim <- filter(ovi_roh_go, Benjamini <= 0.05)
botat_roh_go_slim <- filter(botat_roh_go, Benjamini <= 0.05)
common_go_slim <- filter(common_go, Benjamini <= 0.05)
ev_a_go_slim <- filter(ev_a_go, Benjamini <= 0.05)
ev_b_go_slim <- filter(ev_b_go, Benjamini <= 0.05)
ev_c_go_slim <- filter(ev_c_go, Benjamini <= 0.05)
ovi_go_slim <- filter(ovi_go, Benjamini <= 0.05)
botat_go_slim <- filter(botat_go, Benjamini <= 0.05)

ovi_cnv_go_plot <- ggplot(ovi_cnv_go_slim, aes(x=Fold.enrichment, y=reorder(ID, -Fold.enrichment), fill = Benjamini, size = Result.count)) +
  geom_point(shape=21) +
  scale_fill_viridis(discrete=F, option="magma") +
  scale_size(range = c(3, 5), name = "Result count") +
  theme_minimal() +
  ggtitle ("OVI CNV") +
  ylab("") +
  xlab("Fold enrichment") +
  theme(text = element_text(size = 15)) 

ovi_roh_go_plot <- ggplot(ovi_roh_go_slim, aes(x=Fold.enrichment, y=reorder(ID, -Fold.enrichment), fill = Benjamini, size = Result.count)) +
  geom_point(shape=21) +
  scale_fill_viridis(discrete=F, option="magma") +
  scale_size(range = c(3, 5), name = "Result count") +
  theme_minimal() +
  ggtitle ("OVI ROH") +
  ylab("") +
  xlab("Fold enrichment") +
  theme(text = element_text(size = 15)) 

botat_roh_go_plot <- ggplot(botat_roh_go_slim, aes(x=Fold.enrichment, y=reorder(ID, -Fold.enrichment), fill = Benjamini, size = Result.count)) +
  geom_point(shape=21) +
  scale_fill_viridis(discrete=F, option="magma") +
  scale_size(range = c(3, 5), name = "Result count") +
  theme_minimal() +
  ggtitle ("BoTat ROH") +
  ylab("") +
  xlab("Fold enrichment") +
  theme(text = element_text(size = 15)) 

common_go_plot <- ggplot(common_go_slim, aes(x=Fold.enrichment, y=reorder(ID, -Fold.enrichment), fill = Benjamini, size = Result.count)) +
  geom_point(shape=21) +
  scale_fill_viridis(discrete=F, option="magma") +
  scale_size(range = c(2, 6), name = "Result count") +
  theme_minimal() +
  ggtitle ("Common dN/dS") +
  ylab("") +
  xlab("Fold enrichment") +
  theme(text = element_text(size = 15)) 

ev_a_go_plot <- ggplot(ev_a_go_slim, aes(x=Fold.enrichment, y=reorder(ID, -Fold.enrichment), fill = Benjamini, size = Result.count)) +
  geom_point(shape=21) +
  scale_fill_viridis(discrete=F, option="magma") +
  scale_size(range = c(1,5), name = "Result count") +
  theme_minimal() + 
  ggtitle ("ev.A dN/dS") +
  ylab("") +
  xlab("Fold enrichment") +
  theme(text = element_text(size = 15)) 

ev_b_go_plot <- ggplot(ev_b_go_slim, aes(x=Fold.enrichment, y=reorder(ID, -Fold.enrichment), fill = Benjamini, size = Result.count)) +
  geom_point(shape=21) +
  scale_fill_viridis(discrete=F, option="magma") +
  scale_size(range = c(2, 6), name = "Result count") +
  theme_minimal() + 
  ggtitle ("ev. B dN/dS") +
  ylab("") +
  xlab("Fold enrichment") +
  theme(text = element_text(size = 15)) 

ev_c_go_plot <- ggplot(ev_c_go_slim, aes(x=Fold.enrichment, y=reorder(ID, -Fold.enrichment), fill = Benjamini, size = Result.count)) +
  geom_point(shape=21) +
  scale_fill_viridis(discrete=F, option="magma") +
  scale_size(range = c(2, 6), name = "Result count") +
  theme_minimal() +
  ggtitle ("ev. IVM-t1 dN/dS") +
  ylab("") +
  xlab("Fold enrichment") +
  theme(text = element_text(size = 15)) 

ovi_go_plot <- ggplot(ovi_go_slim, aes(x=Fold.enrichment, y=reorder(ID, -Fold.enrichment), fill = Benjamini, size = Result.count)) +
  geom_point(shape=21) +
  scale_fill_viridis(discrete=F, option="magma") +
  scale_size(range = c(2, 6), name = "Result count") +
  theme_minimal() + 
  ggtitle ("OVI dN/dS") +
  ylab("") +
  xlab("Fold enrichment") +
  theme(text = element_text(size = 15)) 

botat_go_plot <- ggplot(botat_go_slim, aes(x=Fold.enrichment, y=reorder(ID, -Fold.enrichment), fill = Benjamini, size = Result.count)) +
  geom_point(shape=21) +
  scale_fill_viridis(discrete=F, option="magma") +
  scale_size(range = c(2, 6), name = "Result count") +
  theme_minimal() + 
  ggtitle ("BoTat dN/dS") +
  ylab("") +
  xlab("Fold enrichment") +
  theme(text = element_text(size = 15)) 

#Supplementary figure S1
setwd("/Volumes/matthews/Guy/Raw_data/monomorph/data/codeml")
dnds <- read.csv("dnds.txt",header=TRUE)

keycol <- "Gene ID"
valuecol <- "dnds"
gathercols <- c("T.b.evansi.type.A", "T.b.evansi.type.B", "T.b.evansi.type.IVM.t1", "T.b.equiperdum.type.OVI", "T.b.equiperdum.type.BoTat")

dnds_long <- gather(dnds, keycol, valuecol, gathercols)
dnds_long$keycol <- sub('T.b.evansi.type.A','ev.A', dnds_long$keycol)
dnds_long$keycol <- sub('T.b.evansi.type.B','ev.B', dnds_long$keycol)
dnds_long$keycol <- sub('T.b.evansi.type.IVM.t1','ev.IVM-t1', dnds_long$keycol)
dnds_long$keycol <- sub('T.b.equiperdum.type.OVI','OVI',dnds_long$keycol)
dnds_long$keycol <- sub('T.b.equiperdum.type.BoTat','BoTat',dnds_long$keycol)

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}
dnds_long$change <- dnds_long$valuecol-dnds_long$background

dnds_filtered <- dnds_long %>%  filter(change > -50)

dnds_plot <- ggplot(dnds_long, aes(x=(valuecol), y=keycol, colour=Category)) +
  geom_density_ridges(alpha = 0) +
  scale_y_discrete(expand = expansion(add = c(.3, 1.8))) +
  theme_bw() + 
  ggtitle ("dN/dS") +
  ylab("") +
  xlab("dN/dS") +
  theme(text = element_text(size = 15)) 

tiff("/Users/goldriev/Google Drive/My Drive/Developmental_competence_ms/draft_ms/figures/Fig.1/Fig.S1.tiff", units="in", width=15, height=20, res=300)
ggarrange(ovi_cnv_go_plot, ev_a_go_plot, ev_b_go_plot, ev_c_go_plot, ovi_go_plot, botat_go_plot, dnds_plot, ovi_roh_go_plot, ncol = 2, nrow = 4, common.legend = F, legend="right", align = c("h"), labels = "auto", font.label = list(size = 20, color = "black", face = "bold", family = NULL))
dev.off()
