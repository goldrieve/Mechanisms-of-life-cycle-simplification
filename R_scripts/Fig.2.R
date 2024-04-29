library(ggpubr)
library(RColorBrewer)
library(rstatix)
library(viridis)
library(dplyr) 
library(tidyverse)
library(ggpmisc)
library(reshape2)

#Pre-process data
setwd("/Users/goldriev/keep/Store/PHD/TERGO/mono_variant_calling/CRISPR/BHI")
BHI <- read.csv("mono_replace_GC.csv")
IFA <- read.csv("IF.csv")

names(IFA)[names(IFA) == "X1K1N"] <- "1K1N"
names(IFA)[names(IFA) == "X2K1N"] <- "2K1N"
names(IFA)[names(IFA) == "X2K2N"] <- "2K2N"

IFA$"2K1N/2K2N" <- (((IFA$"2K2N" + IFA$"2K1N")/ IFA$Total)*100)

IFA_out <- split( IFA , f = IFA$Cell.line)

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

pal = c("black", "#E41A1C", "#387EB8")

#Tb927.2.4020

out <- split( BHI , f = BHI$Cell_line )
Tb927.2.4020_df <- rbind(out$`IVM-t1_3_Tb927.2.4020`, out$`WT-2.4020_Tb927.2.4020`, out$`APPBP1_AB_P+G_2_Tb927.2.4020`)
Tb927.2.4020_df$Oligopeptides <- factor(Tb927.2.4020_df$Oligopeptides, levels=c("HMI-9", "HMI-9:BHI"))
Tb927.2.4020_df$clade <- factor(Tb927.2.4020_df$Clade.1, levels=c("brucei", "ev.IVM-t1", "add-back"))

Tb927.2.4020_df %>%
  group_by(clade, Hours) %>%
  get_summary_stats(Norm_density, type = "mean_sd")

Tb927.2.4020_df %>%
  group_by(clade, Hours, Oligopeptides) %>%
  identify_outliers(Norm_density)

res.aov <- Tb927.2.4020_df %>% 
  group_by(Oligopeptides) %>%
  anova_test(dv = Norm_density, wid = flask, within = Hours, between = clade)

pwc <- Tb927.2.4020_df %>%
  group_by (Oligopeptides, Hours) %>%
  pairwise_t_test(
    Norm_density ~ clade, paired = F,
    p.adjust.method = "bonferroni"
  )

pwc <- pwc %>% add_xy_position(x = "Hours")
pwc$test <- paste (pwc$group1, c("vs"), pwc$group2)
pwc_slim <-filter(pwc, test != 'ev.IVM-t1 vs add-back')

a <- ggline(Tb927.2.4020_df, x ="Hours", y = "Norm_density", 
            add = c("mean_se", "jitter"),  
            size=1, 
            add.params = list(size = 2),
            color = "clade", facet.by = "Oligopeptides", palette = pal) +
  scale_y_continuous(label=scientific_10) +
  ylab ("Parasites / ml") + 
  ggtitle ("Tb927.2.4020") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 12)) +  
  stat_pvalue_manual(x = "Hours", pwc_slim,  label = "p.adj.signif", tip.length = 0, hide.ns = T, color = 'group2', size = 5) 

IFA_b <- rbind(IFA_out$`T. b. brucei_Tb927.2.4020`, IFA_out$`T.b.evansi type IVM-t1_Tb927.2.4020`, IFA_out$`T. b. brucei A/B 2_Tb927.2.4020`)
IFA_b  <- melt(IFA_b[,c('spp', 'Line', 'BHI', '2K1N/2K2N','PAD')],)

positions <- c("brucei", "ev.IVM-t1", "add-back")

b <- ggplot(IFA_b, aes(x = spp, y = value)) + 
  geom_bar(aes(fill = variable), colour = "black", stat = "identity", position = "dodge") +
  scale_fill_manual(values=c("#999999", "#E69F00")) +
  ylim(0, 80) +
  facet_wrap(c("BHI", "variable")) +
  theme_minimal() + 
  theme(text = element_text(size = 12)) +
  ylab("%") +
  xlab("") +
  ggtitle ("Tb927.2.4020") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(limits = positions)

#Tb927.5.2580

Tb927.5.2580_df <- rbind(out$WT_1_Tb927.5.2580, out$MU09_1_Tb927.5.2580, out$`5_AB_G418/Phleo_6_Tb927.5.2580`)
Tb927.5.2580_df$Oligopeptides <- factor(Tb927.5.2580_df$Oligopeptides, levels=c("HMI-9", "HMI-9:BHI"))
Tb927.5.2580_df$clade <- factor(Tb927.5.2580_df$Clade.1, levels=c("brucei", "ev.A", "add-back"))

Tb927.5.2580_df %>%
  group_by(clade, Hours) %>%
  get_summary_stats(Norm_density, type = "mean_sd")

Tb927.5.2580_df %>%
  group_by(clade, Hours, Oligopeptides) %>%
  identify_outliers(Norm_density)

res.aov <- Tb927.5.2580_df %>% 
  group_by(Oligopeptides) %>%
  anova_test(dv = Norm_density, wid = flask, within = Hours, between = clade)

pwc <- Tb927.5.2580_df %>%
  group_by (Oligopeptides, Hours) %>%
  pairwise_t_test(
    Norm_density ~ clade, paired = F,
    p.adjust.method = "bonferroni"
  )

pwc <- pwc %>% add_xy_position(x = "Hours")
pwc$test <- paste (pwc$group1, c("vs"), pwc$group2)
pwc_slim <-filter(pwc, test != 'ev.A vs add-back')

c <- ggline(Tb927.5.2580_df, x ="Hours", y = "Norm_density", 
       add = c("mean_se", "jitter"),  
       size=1, 
       add.params = list(size = 2),
       color = "clade", facet.by = "Oligopeptides", palette = pal) +
  scale_y_continuous(label=scientific_10) +
  ylab ("Parasites / ml") + 
  ggtitle ("Tb927.5.2580") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 12)) +
  stat_pvalue_manual(x = "Hours", pwc_slim,  label = "p.adj.signif", tip.length = 0, hide.ns = T, color = 'group2', size = 5) 

IFA_d <- rbind(IFA_out$`T. b. brucei_Tb927.5.2580`, IFA_out$`T. b. evansi type A_Tb927.5.2580`, IFA_out$`Add-back evansi A/ evansi A_Tb927.5.2580`)

IFA_d  <- melt(IFA_d[,c('spp', 'Line', 'BHI', '2K1N/2K2N','PAD')],)

positions <- c("brucei", "ev.A", "add-back")

d <- ggplot(IFA_d, aes(x = spp, y = value)) + 
  geom_bar(aes(fill = variable), colour = "black", stat = "identity", position = "dodge") +
  scale_fill_manual(values=c("#999999", "#E69F00")) +
  ylim(0, 80) +
  facet_wrap(c("BHI", "variable")) +
  theme_minimal() + 
  theme(text = element_text(size = 12)) +
  ylab("%") +
  xlab("") +
  ggtitle ("Tb927.5.2580") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "right") +
  scale_x_discrete(limits = positions)

Tb927.11.3400_df <- rbind(out$'11_WT_3_Tb927.11.3400', out$'11_OVI_2_Tb927.11.3400', out$'11_BoTat_2_Tb927.11.3400')
Tb927.11.3400_df$Oligopeptides <- factor(Tb927.11.3400_df$Oligopeptides, levels=c("HMI-9", "HMI-9:BHI"))
Tb927.11.3400_df$clade <- factor(Tb927.11.3400_df$Clade.1, levels=c("brucei", "BoTat", "OVI"))

Tb927.11.3400_df %>%
  group_by(clade, Hours) %>%
  get_summary_stats(Norm_density, type = "mean_sd")

Tb927.11.3400_df %>%
  group_by(clade, Hours, Oligopeptides) %>%
  identify_outliers(Norm_density)

res.aov <- Tb927.11.3400_df %>% 
  group_by(Oligopeptides) %>%
  anova_test(dv = Norm_density, wid = flask, within = Hours, between = clade)

pwc <- Tb927.11.3400_df %>%
  group_by (Oligopeptides, Hours) %>%
  pairwise_t_test(
    Norm_density ~ clade, paired = F,
    p.adjust.method = "bonferroni"
  )

pwc <- pwc %>% add_xy_position(x = "Hours")
pwc$test <- paste (pwc$group1, c("vs"), pwc$group2)
pwc_slim <-filter(pwc, test != 'BoTat vs OVI')

e <- ggline(Tb927.11.3400_df, x ="Hours", y = "Norm_density", 
       add = c("mean_se", "jitter"),  
       size=1, 
       add.params = list(size = 2),
       color = "clade", facet.by = "Oligopeptides", palette = pal) +
  scale_y_continuous(label=scientific_10) +
  ylab ("Parasites / ml") + 
  ggtitle ("Tb927.11.3400") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 12)) +
  stat_pvalue_manual(x = "Hours", pwc_slim,  label = "p.adj.signif", tip.length = 0, hide.ns = T, color = 'group2', size = 5) 

#Motility plot
setwd("/Users/goldriev/keep/Store/PHD/TERGO/mono_variant_calling/CRISPR/FAZ41_motility/")
motility <- read.csv("mash.csv")
motility$SPP <- factor(motility$SPP, levels=c("brucei", "BoTat", "OVI", "BoTat:add-back", "OVI:add-back"))

faz_comparisons <- list( c("brucei", "OVI"), 
                         c("brucei", "BoTat"),  
                         c("brucei", "BoTat:add-back"),  
                         c("brucei", "OVI:add-back"),
                         c("BoTat", "BoTat:add-back"),  
                         c("OVI", "OVI:add-back")  )

pal = c("black", "#E41A1C", "#387EB8","#29AF7FFF", "pink")

f <- ggboxplot(motility, x = "SPP", y = "MEAN_STRAIGHT_LINE_SPEED",
                 add.params = list(size = 2, alpha = 0.1), color = "SPP", palette = pal,
                 add = "jitter") +
  stat_compare_means(comparisons = faz_comparisons, label = ("p.signif")) +
  ylab ("Mean straight line speed (microns/s)") +
  xlab ("Clade") +
  ggtitle ("Tb927.11.3400") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 12))  + 
  guides(color = FALSE, size = FALSE) +
  labs(fill="Clade")

tiff("/Users/s1886853/Google Drive/My Drive/Developmental_competence_ms/draft_ms/figures/Fig.2/Fig.2.tiff", units="in", width=13, height=15, res=300)
ggarrange(a, b, c, d, e, f, ncol = 2, nrow = 3, common.legend = F, legend="top", align = c("hv"), labels = "auto", font.label = list(size = 14, color = "black", face = "bold", family = NULL, padding = unit(1000,"line"))) 
dev.off()

#Supplementary

#Tb927.2.4020

Tb927.2.4020_supp_df <- rbind(out$BoTat_Tb927.2.4020, out$MU09_Tb927.2.4020, out$`WT-2.4020_Tb927.2.4020`)
Tb927.2.4020_supp_df$Oligopeptides <- factor(Tb927.2.4020_supp_df$Oligopeptides, levels=c("HMI-9", "HMI-9:BHI"))

sa <- ggline(Tb927.2.4020_supp_df, x ="Hours", y = "Norm_density", 
             add = c("mean_se", "jitter"),  
             size=1, 
             add.params = list(size = 2),
             color = "Clade.1", facet.by = "Oligopeptides", palette = pal) +
  ggtitle ("Tb927.2.4020") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(label=scientific_10) +
  ylab ("Parasites / ml") + 
  labs(color="Clade") +
  theme(text = element_text(size = 12)) 

#Tb927.8.1530

Tb927.8.1530_df <- rbind(out$'8_BoTat_1_Tb927.8.1530', out$'8_WT_5_Tb927.8.1530')
Tb927.8.1530_df$Oligopeptides <- factor(Tb927.8.1530_df$Oligopeptides, levels=c("HMI-9", "HMI-9:BHI"))

sb <- ggline(Tb927.8.1530_df, x ="Hours", y = "Norm_density", 
             add = c("mean_se", "jitter"),  
             size=1, 
             add.params = list(size = 2),
             color = "Clade.1", facet.by = "Oligopeptides", palette = pal) +
  ggtitle ("Tb927.8.1530") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(label=scientific_10) +
  ylab ("Parasites / ml") + 
  labs(color="Clade") +
  theme(text = element_text(size = 12)) 

#Tb927.11.6600

Tb927.11.6600_df <- rbind(out$"11_BoTat_2_Tb927.11.3400", out$"11_OVI_2_Tb927.11.3400", out$"11_WT_3_Tb927.11.3400")
Tb927.11.6600_df$Oligopeptides <- factor(Tb927.11.6600_df$Oligopeptides, levels=c("HMI-9", "HMI-9:BHI"))

sc <- ggline(Tb927.11.6600_df, x ="Hours", y = "Norm_density", 
             add = c("mean_se", "jitter"),  
             size=1, 
             add.params = list(size = 2),
             color = "Clade.1", facet.by = "Oligopeptides", palette = pal) +
  ggtitle ("Tb927.11.6600") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(label=scientific_10) +
  ylab ("Parasites / ml") + 
  labs(color="Clade") +
  theme(text = element_text(size = 12)) 


#Tb927.4.3650 

Tb927.4.3650_df <- rbind(out$"4_WT_1_Tb927.4.3560", out$"4_MU10_2_Tb927.4.3560")
Tb927.4.3650_df$Oligopeptides <- factor(Tb927.4.3650_df$Oligopeptides, levels=c("HMI-9", "HMI-9:BHI"))

sd <- ggline(Tb927.4.3650_df, x ="Hours", y = "Norm_density", 
             add = c("mean_se", "jitter"),  
             size=1, 
             add.params = list(size = 2),
             color = "Clade.1", facet.by = "Oligopeptides", palette = pal) +
  ggtitle ("Tb927.4.3650") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(label=scientific_10) +
  ylab ("Parasites / ml") + 
  labs(color="Clade") +
  theme(text = element_text(size = 12)) 

#G224S

G224S_df <- rbind(out$`G224S+/+_Tb927.2.4020`, out$`WT-2.4020_Tb927.2.4020`)
G224S_df$Oligopeptides <- factor(G224S_df$Oligopeptides, levels=c("HMI-9", "HMI-9:BHI"))
G224S_df$Clade <- factor(G224S_df$Clade, levels=c("T.b.brucei", "G224S"))

se <- ggline(G224S_df, x ="Hours", y = "Norm_density", 
             add = c("mean_se", "jitter"),  
             size=1, 
             add.params = list(size = 2),
             color = "Clade.1", facet.by = "Oligopeptides", palette = pal) +
  ggtitle ("G224S") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(label=scientific_10) +
  ylab ("Parasites / ml") + 
  labs(color="Clade") +
  theme(text = element_text(size = 12)) 

#A149P

A149P_df <- rbind(out$`A149P_H/B_5_Tb927.5.2580`, out$WT_1_Tb927.5.2580)
A149P_df$Oligopeptides <- factor(A149P_df$Oligopeptides, levels=c("HMI-9", "HMI-9:BHI"))
A149P_df$Clade <- factor(A149P_df$Clade, levels=c("T.b.brucei", "A149P"))

sf <- ggline(A149P_df, x ="Hours", y = "Norm_density", 
       add = c("mean_se", "jitter"),  
       size=1, 
       add.params = list(size = 2),
       color = "Clade.1", facet.by = "Oligopeptides", palette = pal) +
  ggtitle ("A149P") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(label=scientific_10) +
  ylab ("Parasites / ml") + 
  labs(color="Clade") +
  theme(text = element_text(size = 12)) 

tiff("/Users/s1886853/Google Drive/My Drive/Developmental_competence_ms/draft_ms/figures/Fig.2/S3.tiff", units="in", width=12, height=12, res=300)
ggarrange(sa, sb, sc, sd, se, sf, ncol = 2, nrow = 3, common.legend = F, legend="top", align = c("hv"), labels = "auto", font.label = list(size = 14, color = "black", face = "bold", family = NULL))
dev.off()
