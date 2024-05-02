# Define libraries
libraries <- c("ggplot2", "ggrepel", "RColorBrewer", "ggpubr", "dplyr", "ggridges", 
               "splitstackshape", "tidyr", "circlize", "viridis", 
               "ggplotify", "gridBase", "ape", "Biostrings", "ggtree", "treeio", 
               "phytools", "cowplot", "gridExtra")
# Load   
invisible(lapply(libraries, library, character.only = TRUE))

# Set up working directory
setwd("/Volumes/matthews/Guy/Raw_data/monomorph/data/circos")

# Read in chromosome data
chr <- read.table('chr.txt')

# Define files
file_lists <- list(
  shuf_files = c('T.b.evansitypeA.vcf.variants.shuf', 'T.b.evansitypeB.vcf.variants.shuf', 
                 'T.b.evansitypeC.vcf.variants.shuf', 'T.b.equitypeOVI.vcf.variants.shuf', 
                 'T.b.equitypeBOTAT.vcf.variants.shuf'),
  roh_files = c('evansi_a_roh.txt', 'evansi_b_roh.txt', 
                'evansi_c_roh.txt', 'equi_ovi_roh.txt', 
                'equi_botat_roh.txt'),
  cnv_files = c('evansi_a_cnv.txt', 'evansi_b_cnv.txt', 
                'evansi_c_cnv.txt', 'equi_ovi_cnv.txt', 
                'equi_botat_cnv.txt'),
  dnds_files = c('evansi_a_dnds.txt', 'evansi_b_dnds.txt', 
                 'evansi_c_dnds.txt', 'equi_ovi_dnds.txt', 
                 'equi_botat_dnds.txt')
)

# Function to read files
read_files <- function(files) {
  lapply(files, read.table, header = F)
}

# Read files
file_data <- lapply(file_lists, read_files)

# Plot circos 
circlize_plot = function() {
  circos.par("track.height" = 0.15, cell.padding = c(0, 0, 0, 0), "clock.wise" = T, "gap.after" = 10)
  circos.par(gap.after = c(rep(1, 10), 20))
  circos.genomicInitialize(chr,sector.names = c("1","2", "3", "4", "5", "6", "7", "8", "9", "10", "11"), labels.cex = 0.8)
  col_mat <-  c("darkred", "darkgreen", "#643F95", "gold", "darkorange")
  circos.genomicDensity(file_data$shuf_files, col = col_mat, baseline = 0, area = F, window.size = 200000, count_by = c("percent"))
  f = colorRamp2(breaks = c(1.9, 2.1), colors = c("darkblue", "red"))
  
  for (x in file_data$cnv_files) {
    circos.genomicTrack(x, ylim =c(0,0.01), track.height = 0.02, panel.fun = function(region, value, ...) {
      i = getI(...)
      circos.genomicRect(region, value, col = f(value[[1]]), border = NA)})}
  
  circos.genomicDensity(file_data$dnds_files, col = col_mat, baseline = 0, area = F, window.size = 200000, count_by = c("percent"))
  
  for (x in file_data$roh_files) {
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

# Add text to circos plot
text(.27, .06, substitute(paste(bold("ROH"))), family = "helvetica", cex = 0.6)
text(.43, .06, substitute(paste(bold("dN/dS"))), family = "helvetica", cex = 0.6)
text(.63, .06, substitute(paste(bold("CNV"))), family = "helvetica", cex = 0.6)
text(.81, .06, substitute(paste(bold("SNP"))), family = "helvetica", cex = 0.6)

# Draw image to create ggplot
p1 <- recordPlot() 
ggdraw(p1)

# Load tree file 
setwd("/Volumes/matthews/Guy/Raw_data/monomorph/data/figures/Fig.1")
tree <- treeio::read.newick("tree.txt", node.label="support")
root <- rootnode(tree)  

# Plot tree
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

# Export plot
tiff("/Users/goldriev/Google Drive/My Drive/Developmental_competence_ms/draft_ms/figures/Fig.1/Fig.1.tiff", units="in", width=6, height=10, res=300)
ggarrange(plotlist = fig1_plots, ncol = 1, nrow = 2, labels = c("b", "a"), font.label = list(size = 14, color = "black", face = "bold", family = NULL))
dev.off()

# Load gene ontology data for S1

setwd("/Volumes/matthews/Guy/Raw_data/monomorph/data/go_term/")

# Function to read, filter and plot data
process_data <- function(file_name, title) {
  data <- read.csv(file_name)
  data_slim <- filter(data, Benjamini <= 0.05)
  
  plot <- ggplot(data_slim, aes(x=Fold.enrichment, y=reorder(ID, -Fold.enrichment), fill = Benjamini, size = Result.count)) +
    geom_point(shape=21) +
    scale_fill_viridis(discrete=F, option="magma") +
    scale_size(range = c(2, 6), name = "Result count") +
    theme_minimal() +
    ggtitle (title) +
    ylab("") +
    xlab("Fold enrichment") +
    theme(text = element_text(size = 15))
  
  return(plot)
}

# Set working directory
setwd("/Volumes/matthews/Guy/Raw_data/monomorph/data/go_term/")

# List of file names and corresponding titles
file_names <- c("ovi_cnv_go.csv", "botat_roh_go.csv", "ovi_roh_go.csv", "common_dnds_go.csv", "evansi_a_dnds_go.csv", "evansi_b_dnds_go.csv", "evansi_c_dnds_go.csv", "equi_ovi_dnds_go.csv", "equi_botat_dnds_go.csv")
titles <- c("OVI CNV", "BoTat ROH", "OVI ROH", "Common dN/dS", "ev.A dN/dS", "ev. B dN/dS", "ev. IVM-t1 dN/dS", "OVI dN/dS", "BoTat dN/dS")

# Plot for each file
plots <- mapply(process_data, file_name = file_names, title = titles, SIMPLIFY = FALSE)

#Supplementary figure S1
setwd("/Volumes/matthews/Guy/Raw_data/monomorph/data/codeml")
dnds <- read.csv("dnds.txt",header=TRUE)

# Set up data structure
keycol <- "Gene ID"
valuecol <- "dnds"
gathercols <- c("T.b.evansi.type.A", "T.b.evansi.type.B", "T.b.evansi.type.IVM.t1", "T.b.equiperdum.type.OVI", "T.b.equiperdum.type.BoTat")

# Make a long form dataframe
dnds_long <- dnds %>%
  gather(keycol, valuecol, gathercols) %>%
  mutate(keycol = sub('T.b.evansi.type.A', 'ev.A', 
                      sub('T.b.evansi.type.B', 'ev.B', 
                          sub('T.b.evansi.type.IVM.t1', 'ev.IVM-t1', 
                              sub('T.b.equiperdum.type.OVI', 'OVI',
                                  sub('T.b.equiperdum.type.BoTat', 'BoTat', keycol))))))

# Plot dnds
dnds_plot <- ggplot(dnds_long, aes(x=(valuecol), y=keycol, colour=Category)) +
  geom_density_ridges(alpha = 0) +
  scale_y_discrete(expand = expansion(add = c(.3, 1.8))) +
  theme_bw() + 
  ggtitle ("dN/dS") +
  ylab("") +
  xlab("dN/dS") +
  theme(text = element_text(size = 15)) 

# Export plot
tiff("/Users/goldriev/Google Drive/My Drive/Developmental_competence_ms/draft_ms/figures/Fig.1/Fig.S1.tiff", units="in", width=15, height=20, res=300)
ggarrange(plots$ovi_cnv_go.csv, plots$evansi_a_dnds_go.csv, plots$evansi_b_dnds_go.csv, plots$evansi_c_dnds_go.csv, plots$equi_ovi_dnds_go.csv, plots$equi_botat_dnds_go.csv, dnds_plot, plots$ovi_roh_go.csv, ncol = 2, nrow = 4, common.legend = F, legend="right", align = c("h"), labels = "auto", font.label = list(size = 20, color = "black", face = "bold", family = NULL))
dev.off()