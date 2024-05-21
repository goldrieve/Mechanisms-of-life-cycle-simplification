# Libraries
libraries <- c("ggplot2", "ggpubr", "ggrepel", "tidyverse", "rstatix", "dplyr", "viridis", "gridExtra", "ggpp", "ggpmisc")

# Load   
invisible(lapply(libraries, library, character.only = TRUE))

# Load data and calculate stats
setwd("/Volumes/matthews/Guy/Raw_data/monomorph/data/figures/Fig.3")
vivo <- read.csv("Tb927.5.2580.csv")
vivo$PAD_P <- (vivo$PAD/(vivo$X1K1N + vivo$X2K1N + vivo$X2K2N + vivo$Other))*100
vivo$dividing <- ((vivo$X2K1N + vivo$X2K2N)/ vivo$Total)*100
vivo$Clade <- factor(vivo$Clade, levels=c("brucei", "ev.A", "add-back"))

# Set up
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

pal = c("black", "#E41A1C", "#387EB8")

# Analyse 
vivo %>%
  group_by(Clade, DPI) %>%
  get_summary_stats(cells.ml, type = "mean_sd")

vivo %>%
  group_by(Clade, DPI) %>%
  identify_outliers(cells.ml)

res.aov <- vivo %>% 
  anova_test(dv = cells.ml, wid = ID, within = DPI, between = Clade)

pwc <- vivo %>%
  group_by(DPI) %>%
  pairwise_t_test(
    cells.ml ~ Clade, paired = F,
    p.adjust.method = "bonferroni"
  )

pwc <- pwc %>% add_xy_position(x = "DPI")
pwc$test <- paste (pwc$group1, c("vs"), pwc$group2)
pwc_slim <- filter(pwc, test != 'ev.A vs add-back')
pwc_slim$xmin <- as.numeric(c("4", "4", "5", "5", "6", "6", "7", "7", "8", "8", "9", "9"))
pwc_slim$xmax <- as.numeric(c("4", "4", "5", "5", "6", "6", "7", "7", "8", "8", "9", "9"))

cells <- ggline(vivo, x = "DPI", y = "cells.ml", 
                add = c("mean_se", "jitter"),  
                size=1, 
                add.params = list(size = 2),
                color = "Clade", palette = pal) +
  ggtitle("Tb927.5.2580 parasitemia") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab ("Parasites / ml") + 
  theme(text = element_text(size = 15)) +
  scale_y_continuous(label=scientific_10) +
  stat_pvalue_manual(pwc_slim, label = "p.adj.signif", tip.length = 0, hide.ns = T, color = 'group2', size = 5) 

cells_mouse <- ggline(vivo, x = "DPI", y = "cells.ml",  
       size=1, 
       add.params = list(size = 2.5),
       color = "Clade", 
       facet.by = "Mouse", 
       short.panel.labs = T, palette = pal) +
  #scale_color_viridis(discrete = T) +
  ggtitle("Tb927.5.2580 parasitemia") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab ("DPI") +
  ylab ("Parasites / ml") + 
  theme(text = element_text(size = 15)) +
  scale_y_continuous(label=scientific_10) 

#PAD 

filtered <- vivo %>% 
  filter(DPI!="4")

filtered %>%
  group_by(Clade, DPI) %>%
  get_summary_stats(PAD_P, type = "mean_sd")

filtered %>%
  group_by(Clade, DPI) %>%
  identify_outliers(PAD_P)

res.aov <- filtered %>% 
  anova_test(dv = PAD_P, wid = ID, within = DPI, between = Clade)

pwc <- filtered %>%
  group_by(DPI) %>%
  pairwise_t_test(
    PAD_P ~ Clade, paired = F,
    p.adjust.method = "bonferroni"
  )

pwc <- pwc %>% add_xy_position(x = "DPI")
pwc$test <- paste (pwc$group1, c("vs"), pwc$group2)
pwc_slim <-filter(pwc, test != 'ev.A vs add-back')
pwc_slim$xmin <- as.numeric(c("5", "5", "6", "6", "7", "7", "8", "8", "9", "9"))
pwc_slim$xmax <- as.numeric(c("5", "5", "6", "6", "7", "7", "8", "8", "9", "9"))

pad <- ggline(filtered, x = "DPI", y = "PAD_P", 
               add = c("mean_se", "jitter"),  
               size=1, 
               add.params = list(size = 2),
               color = "Clade", palette = pal) +
  ggtitle("Tb927.5.2580 PAD1+") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab ("PAD1+ (%)") + 
  xlim (4,9) +
  theme(text = element_text(size = 15)) +
  stat_pvalue_manual(pwc_slim,  label = "p.adj.signif", tip.length = 0, hide.ns = T, color = 'group2', size = 5) 

pad_mouse <- ggline(vivo, x = "DPI", y = "PAD_P",  
                   size=1, 
                   add.params = list(size = 2.5),
                   color = "Clade", 
                   facet.by = "Mouse", 
                   short.panel.labs = T, palette = pal) +
  #scale_color_viridis(discrete = T) +
  ggtitle("Tb927.5.2580 PAD1+") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab ("DPI") +
  ylab ("PAD1+ (%)") + 
  xlim (4,9) +
  theme(text = element_text(size = 15)) 

#dividing

filtered %>%
  group_by(Clade, DPI) %>%
  get_summary_stats(dividing, type = "mean_sd")

filtered %>%
  group_by(Clade, DPI) %>%
  identify_outliers(dividing)

res.aov <- filtered %>% 
  anova_test(dv = dividing, wid = ID, within = DPI, between = Clade)

pwc <- filtered %>%
  group_by(DPI) %>%
  pairwise_t_test(
    dividing ~ Clade, paired = F,
    p.adjust.method = "bonferroni"
  )

pwc <- pwc %>% add_xy_position(x = "DPI")
pwc$test <- paste (pwc$group1, c("vs"), pwc$group2)
pwc_slim <-filter(pwc, test != 'ev.A vs add-back')
pwc_slim$xmin <- as.numeric(c("5", "5", "6", "6", "7", "7", "8", "8", "9", "9"))
pwc_slim$xmax <- as.numeric(c("5", "5", "6", "6", "7", "7", "8", "8", "9", "9"))

dividing <- ggline(filtered, x = "DPI", y = "dividing", 
                   add = c("mean_se", "jitter"),  
                   size=1, 
                   add.params = list(size = 3),
                   color = "Clade", palette = pal) +
  #scale_color_viridis(discrete = T) +
  ggtitle("Tb927.5.2580 2K1N/2K2N") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab ("2K1N/2K2N (%)") + 
  theme(text = element_text(size = 15)) +
  xlim (4,9) +
  stat_pvalue_manual(pwc_slim,  label = "p.adj.signif", tip.length = 0, hide.ns = T, color = 'group2', size = 5) 

dividing_mouse <- ggline(vivo, x = "DPI", y = "dividing",  
                         size=1, 
                         add.params = list(size = 2.5),
                         color = "Clade", 
                         facet.by = "Mouse", 
                         short.panel.labs = T, palette = pal) +
  ggtitle("Tb927.5.2580 2K1N/2K2N") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab ("DPI") +
  ylab ("2K1N/2K2N (%)") + 
  xlim (4,9) +
  theme(text = element_text(size = 15)) 

d <- ggplot() + 
  theme_void() + 
  ggtitle("Tb927.5.2580 IFA") + 
  theme(plot.title = element_text(hjust = 0.5))  +
  theme(text = element_text(size = 15)) 


tiff("~/Google Drive/My Drive/Developmental_competence_ms/draft_ms/figures/Fig.3/Fig.3.tiff", units="in", width=12, height=10, res=300)
ggarrange(cells, pad, dividing, d, ncol = 2, nrow = 2, common.legend = F, legend="top", align = c("h"), labels = c("a", "b", "c", "d"), font.label = list(size = 14, color = "black", face = "bold", family = NULL))
dev.off()

tiff("~/Google Drive/My Drive/Developmental_competence_ms/draft_ms/figures/Fig.3/Fig.S5.tiff", units="in", width=12, height=10, res=300)
ggarrange(cells_mouse, pad_mouse, dividing_mouse, ncol = 2, nrow = 2, common.legend = T, legend="top", align = c("h"), labels = "auto", font.label = list(size = 14, color = "black", face = "bold", family = NULL))
dev.off()
