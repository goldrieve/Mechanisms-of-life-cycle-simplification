# Define libraries
libraries <- c("ggpubr", "RColorBrewer", "rstatix", "viridis", "dplyr", "tidyverse", "ggpmisc", "reshape2", "ggh4x")

# Load
invisible(lapply(libraries, library, character.only = TRUE))

# Pre-process data
setwd("/Volumes/matthews/Guy/Raw_data/monomorph/data/figures/Fig.2")
BHI <- read.csv("mono_replace_GC.csv")
IFA <- read.csv("IF.csv")

# Edit column header
name_changes <- c("X1K1N" = "1K1N", "X2K1N" = "2K1N", "X2K2N" = "2K2N")
names(IFA) <- ifelse(names(IFA) %in% names(name_changes), name_changes[names(IFA)], names(IFA))

# Calculate % 2K1N/ 2K2N
IFA$"2K1N/2K2N" <- (((IFA$"2K2N" + IFA$"2K1N")/ IFA$Total)*100)

IFA_out <- split( IFA , f = IFA$Cell.line)

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

pal = c("black", "#E41A1C", "#387EB8")

out <- split( BHI , f = BHI$Cell_line )

# Define a function to perform the operations
analyze_data <- function(df, clade_levels, pwc_filter, plot_title) {
  df$Oligopeptides <- factor(df$Oligopeptides, levels=c("HMI-9", "HMI-9:BHI"))
  df$clade <- factor(df$Clade.1, levels=clade_levels)
  
  df %>%
    group_by(clade, Hours) %>%
    get_summary_stats(Norm_density, type = "mean_sd")
  
  df %>%
    group_by(clade, Hours, Oligopeptides) %>%
    identify_outliers(Norm_density)
  
  res.aov <- df %>% 
    group_by(Oligopeptides) %>%
    anova_test(dv = Norm_density, wid = flask, within = Hours, between = clade)
  
  pwc <- df %>%
    group_by (Oligopeptides, Hours) %>%
    pairwise_t_test(
      Norm_density ~ clade, paired = F,
      p.adjust.method = "bonferroni"
    )
  
  pwc <- pwc %>% add_xy_position(x = "Hours")
  pwc$test <- paste (pwc$group1, c("vs"), pwc$group2)
  pwc_slim <-filter(pwc, test != pwc_filter)
  
  plot <- ggline(df, x ="Hours", y = "Norm_density", 
                 add = c("mean_se", "jitter"),  
                 size=1, 
                 add.params = list(size = 2),
                 color = "clade", facet.by = "Oligopeptides", palette = pal) +
    scale_y_continuous(label=scientific_10) +
    ylab ("Parasites / ml") + 
    ggtitle (plot_title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(face = "bold")) +
    theme(text = element_text(size = 12)) +
    stat_pvalue_manual(x = "Hours", pwc_slim,  label = "p.adj.signif", tip.length = 0, hide.ns = T, color = 'group2', size = 5) 
  
  return(list(res.aov = res.aov, pwc = pwc, plot = plot))
}

# Call the function for each data frame
Tb927.2.4020_df <- rbind(out$`IVM-t1_3_Tb927.2.4020`, out$`WT-2.4020_Tb927.2.4020`, out$`APPBP1_AB_P+G_2_Tb927.2.4020`)
results_2.4020 <- analyze_data(Tb927.2.4020_df, c("brucei", "ev.IVM-t1", "add-back"), 'ev.IVM-t1 vs add-back', "Tb927.2.4020")

Tb927.5.2580_df <- rbind(out$WT_1_Tb927.5.2580, out$MU09_1_Tb927.5.2580, out$`5_AB_G418/Phleo_6_Tb927.5.2580`)
results_5.2580 <- analyze_data(Tb927.5.2580_df, c("brucei", "ev.A", "add-back"), 'ev.A vs add-back', "Tb927.5.2580")

Tb927.11.3400_df <- rbind(out$'11_WT_3_Tb927.11.3400', out$'11_OVI_2_Tb927.11.3400', out$'11_BoTat_2_Tb927.11.3400')
results_11.3400 <- analyze_data(Tb927.11.3400_df, c("brucei", "BoTat", "OVI"), 'BoTat vs OVI', "Tb927.11.3400")

# Plot IFA data for Tb927.2.4020 and Tb927.5.2580

ifa_plot <- function(IFA_data, plot_title, positions) {
  IFA_data <- melt(IFA_data[,c('spp', 'Line', 'BHI', '2K1N/2K2N','PAD1')],)
  
  strip <- strip_themed(background_x = elem_list_rect(fill = c('#E3F2FD', '#E3F2FD', '#DCEDC8','#DCEDC8')))
  
  plot <- ggplot(IFA_data, aes(x = spp, y = value)) + 
    geom_bar(aes(fill = variable), colour = "black", stat = "identity", position = "dodge") +
    scale_fill_manual(values=c("#999999", "#E69F00")) +
    ylim(0, 80) +
    theme_bw() + 
    ylab("%") +
    xlab("") +
    ggtitle (plot_title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(face = "bold")) +
    theme(text = element_text(size = 12)) +
    scale_x_discrete(limits = positions)  +
    facet_wrap2(BHI ~ variable, strip = strip)
  
  return(plot)
}

# Use the function to generate the plots
IFA_b <- rbind(IFA_out$`T. b. brucei_Tb927.2.4020`, IFA_out$`T.b.evansi type IVM-t1_Tb927.2.4020`, IFA_out$`T. b. brucei A/B 2_Tb927.2.4020`)
b <- ifa_plot(IFA_b, "Tb927.2.4020", c("brucei", "ev.IVM-t1", "add-back"))

IFA_d <- rbind(IFA_out$`T. b. brucei_Tb927.5.2580`, IFA_out$`T. b. evansi type A_Tb927.5.2580`, IFA_out$`Add-back evansi A/ evansi A_Tb927.5.2580`)
d <- ifa_plot(IFA_d, "Tb927.5.2580", c("brucei", "ev.A", "add-back"))

#Motility plot

motility <- read.csv("combined_motility.csv")
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
  theme(plot.title = element_text(face = "bold")) +
  theme(text = element_text(size = 12))  + 
  guides(color = FALSE, size = FALSE) +
  labs(fill="Clade") 

tiff("/Users/goldriev/Google Drive/My Drive/Developmental_competence_ms/draft_ms/figures/Fig.2/Fig.2.tiff", units="in", width=13, height=15, res=300)
ggarrange(results_2.4020$plot, b, results_5.2580$plot, d, results_11.3400$plot, f, ncol = 2, nrow = 3, common.legend = F, legend="top", align = c("hv"), labels = "auto", font.label = list(size = 14, color = "black", face = "bold", family = NULL, padding = unit(1000,"line"))) 
dev.off()

#Supplementary

BHI <- read.csv("monomorph_bhi.csv")

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

BHI$change <- ((BHI$Density - BHI$Start)/ BHI$Start)*100
coul <- brewer.pal(5,"Set3") 

s2a <- ggline(BHI, x ="BHI", y = "Density", 
            add = c("mean_se", "jitter"),  
            size=1, 
            add.params = list(size = 2),
            color = "BHI", facet.by = "Line") +
  scale_y_continuous(label = scientific_10) +
  ylab ("Parasites / ml") +
  xlab ("BHI (%)") +
  scale_color_manual(values = coul)


s2b <- ggline(BHI, x ="BHI", y = "change", 
              add = c("mean_se", "jitter"),  
              size=1, 
              add.params = list(size = 2),
              color = "BHI", facet.by = "Line") +
  ylab ("Percentage change") +
  xlab ("BHI (%)") +
  scale_color_manual(values = coul)

s2c <- ggplot() + 
  theme_void() + 
  ggtitle("") 

tiff("/Users/goldriev/Google Drive/My Drive/Developmental_competence_ms/draft_ms/figures/Fig.2/S2.tiff", units="in", width=8, height=14, res=300)
ggarrange(s2a, s2b, s2c, ncol = 1, nrow = 3, common.legend = F, legend="top", align = c("hv"), labels = "auto", font.label = list(size = 14, color = "black", face = "bold", family = NULL, padding = unit(1000,"line"))) 
dev.off()

#Tb927.2.4020

process_data <- function(data_names, title, levels1, levels2 = NULL) {
  df <- do.call(rbind, lapply(data_names, function(x) out[[x]]))
  df$Oligopeptides <- factor(df$Oligopeptides, levels = levels1)
  if (!is.null(levels2)) {
    df$Clade <- factor(df$Clade, levels = levels2)
  }
  
  plot <- ggline(df, x ="Hours", y = "Norm_density", 
                 add = c("mean_se", "jitter"),  
                 size = 1, 
                 add.params = list(size = 2),
                 color = "Clade.1", facet.by = "Oligopeptides", palette = pal) +
    ggtitle (title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(face = "bold")) +
    scale_y_continuous(label = scientific_10) +
    ylab ("Parasites / ml") + 
    labs(color = "Clade") +
    theme(text = element_text(size = 12))
  
  return(plot)
}

sa <- process_data(c("BoTat_Tb927.2.4020", "MU09_Tb927.2.4020", "WT-2.4020_Tb927.2.4020"), "Tb927.2.4020", c("HMI-9", "HMI-9:BHI"))
sb <- process_data(c("8_BoTat_1_Tb927.8.1530", "8_WT_5_Tb927.8.1530"), "Tb927.8.1530", c("HMI-9", "HMI-9:BHI"))
sc <- process_data(c("11_BoTat_2_Tb927.11.3400", "11_OVI_2_Tb927.11.3400", "11_WT_3_Tb927.11.3400"), "Tb927.11.6600", c("HMI-9", "HMI-9:BHI"))
sd <- process_data(c("4_WT_1_Tb927.4.3560", "4_MU10_2_Tb927.4.3560"), "Tb927.4.3650", c("HMI-9", "HMI-9:BHI"))
se <- process_data(c("G224S+/+_Tb927.2.4020", "WT-2.4020_Tb927.2.4020"), "G224S", c("HMI-9", "HMI-9:BHI"), c("T.b.brucei", "G224S"))
sf <- process_data(c("A149P_H/B_5_Tb927.5.2580", "WT_1_Tb927.5.2580"), "A149P", c("HMI-9", "HMI-9:BHI"), c("T.b.brucei", "A149P"))

tiff("/Users/goldriev/Google Drive/My Drive/Developmental_competence_ms/draft_ms/figures/Fig.2/S3.tiff", units="in", width=12, height=12, res=300)
ggarrange(sa, sb, sc, sd, se, sf, ncol = 2, nrow = 3, common.legend = F, legend="top", align = c("hv"), labels = "auto", font.label = list(size = 14, color = "black", face = "bold", family = NULL))
dev.off()