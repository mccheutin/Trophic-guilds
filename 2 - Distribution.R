#### PCOA and nmds on the diet -----
#### 
#### How are distributed the different bacterial composition with the diets?
#### 

setwd("/Users/marie-charlottecheutin/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Carnivory")
dir.create("./Distribution", recursive = T)
setwd("./Distribution")
load("../physeq_object/gut_core.RData")


# Metrics calcul ----
core_rel <- transform_sample_counts(gut_core, function(x) x / sum(x) )
samp_data <- data.frame(sample_data(core_rel))
save(core_rel, file =  "core_rel.RData")

# __ ** Bray-Curtis ----
# PCOA coda & Permanova
library(vegan)
otu.bc <- vegdist(core_rel@otu_table, method = "bray")
pcoa.sub.bc <- pcoa(otu.bc)
pcoa_coord.bc <- pcoa.sub.bc$vectors[,1:3]
# Contruction of the table for graphic 
library(stringr)
hull.bc <- cbind(pcoa_coord.bc, samp_data)
save(hull.bc, file = "hull.bc.RData")
# What is the percentage of the explicative variance? 
paste("Axis 1 :",percent(pcoa.sub.bc$values$Relative_eig[1])) # 10 %
paste("Axis 2 :",percent(pcoa.sub.bc$values$Relative_eig[2])) # 6 %
paste("Axis 3 :",percent(pcoa.sub.bc$values$Relative_eig[3])) # 5 %
paste("Axis 4 :",percent(pcoa.sub.bc$values$Relative_eig[4])) # 4 %
save(otu.bc, hull.bc, file = "beta.bc.RData")

# __ ** Weighted Unifrac ----
# PCOA coda & Permanova
otu.wu <- UniFrac(core_rel, weighted = T, normalized=F, parallel = F, fast=T)
pcoa.sub.wu <- pcoa(otu.wu)
pcoa_coord.wu <- pcoa.sub.wu$vectors[,1:3]
hull.wu <- cbind(pcoa_coord.wu, samp_data)
save(hull.wu, file = "hull.wu.RData")
# What is the percentage of the explicative variance? 
library(scales)
paste("Axis 1 :",percent(pcoa.sub.wu$values$Relative_eig[1])) # 26 %
paste("Axis 2 :",percent(pcoa.sub.wu$values$Relative_eig[2])) # 19 %
paste("Axis 3 :",percent(pcoa.sub.wu$values$Relative_eig[3])) # 12 %
paste("Axis 4 :",percent(pcoa.sub.wu$values$Relative_eig[4])) # 5 %
save(otu.wu, hull.wu, file = "beta.wu.RData")

# __ ** Unweighted Unifrac ----
# PCOA coda & Permanova
otu.uu <- UniFrac(core_rel, weighted = F, normalized=F, parallel = F, fast=T)
pcoa.sub.uu <- pcoa(otu.uu)
pcoa_coord.uu <- pcoa.sub.uu$vectors[,1:3]
hull.uu <- cbind(pcoa_coord.uu, samp_data)
save(hull.uu, file = "hull.uu.RData")

# What is the percentage of the explicative variance? 
library(scales)
paste("Axis 1 :",percent(pcoa.sub.uu$values$Relative_eig[1])) # 14 %
paste("Axis 2 :",percent(pcoa.sub.uu$values$Relative_eig[2])) # 8 %
paste("Axis 3 :",percent(pcoa.sub.uu$values$Relative_eig[3])) # 7 %
paste("Axis 4 :",percent(pcoa.sub.uu$values$Relative_eig[4])) # 5 %
save(otu.uu, hull.uu, file = "beta.uu.RData")

## *** diet_pred *** ----
## __ PCOA plots -----
col_diet_pred <- c('Corallivores' =  'indianred1',
                   'Crustacivores' = 'orange',
                   'Herbivores Microvores Detritivores' = 'darkgreen',
                   'Macroinvertivores' = 'darkslategray',
                   'Microinvertivores' = 'bisque3',
                   'Piscivores' = 'darkblue',
                   'Planktivores' = 'darkred',
                   'sessile invertivores' = 'coral4',
                   'Omnivorous' = 'gray',
                   'Carnivorous' = 'darkslateblue')


pcoa.diet.pred.bc <- ggplot() +
  geom_hline(yintercept=0, colour="lightgrey", linetype = 2) +
  geom_vline(xintercept=0, colour="lightgrey", linetype = 2) +
  scale_color_manual(values = col_diet_pred) +
  geom_point(data = hull.bc, aes(x=Axis.1, y=Axis.2, color = diet_pred), alpha = 0.7, size = 4, shape = 16) +
  xlab(paste("PCo1 (", round(pcoa.sub.bc$values$Relative_eig[1]*100, 1), "%)")) +
  ylab(paste("PCo2 (", round(pcoa.sub.bc$values$Relative_eig[2]*100, 1), "%)"))  +
  theme_bw() +
  #stat_ellipse(data = hull.bc,aes(x=Axis.1, y=Axis.2,color=diet_pred),type = "norm")+
  coord_equal() +
  theme(axis.title.x = element_text(size=18, family = "serif"), # remove x-axis labels
        axis.title.y = element_text(size=18, family = "serif"), # remove y-axis labels
        axis.text.x = element_text(size=14, family = "serif"),
        axis.text.y = element_text(size=14, family = "serif"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.title = element_text(family = "serif", size = 16), 
        legend.position = "none",
        #legend.text = element_text(size = 0, family = "serif"),
        #legend.title = element_text(size = 0,family = "serif"),
        strip.text.x = element_text(size=14, family = "serif"))
#+facet_wrap(~ region)
pcoa.diet.pred.bc

pcoa.diet.pred.wu <- ggplot() +
  geom_hline(yintercept=0, colour="lightgrey", linetype = 2) +
  geom_vline(xintercept=0, colour="lightgrey", linetype = 2) +
  scale_color_manual(values = col_diet_pred) +
  #stat_ellipse(data = hull.wu,aes(x=Axis.1, y=Axis.2,color=diet_pred),type = "norm")+
  geom_point(data = hull.wu, aes(x=Axis.1, y=Axis.2, color = diet_pred), alpha = 0.7, size = 4, shape = 16) +
  xlab(paste("PCo1 (", round(pcoa.sub.wu$values$Relative_eig[1]*100, 1), "%)")) +
  ylab(paste("PCo2 (", round(pcoa.sub.wu$values$Relative_eig[2]*100, 1), "%)"))  +
  theme_bw() +
  coord_equal() +
  theme(axis.title.x = element_text(size=18, family = "serif"), # remove x-axis labels
        axis.title.y = element_text(size=18, family = "serif"), # remove y-axis labels
        axis.text.x = element_text(size=14, family = "serif"),
        axis.text.y = element_text(size=14, family = "serif"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.title = element_text(family = "serif", size = 16), 
        legend.position = "none",
        #legend.text = element_text(size = 0, family = "serif"),
        #legend.title = element_text(size = 0,family = "serif"),
        strip.text.x = element_text(size=14, family = "serif"))
#+facet_wrap(~ region)
pcoa.diet.pred.wu

pcoa.diet.pred.uu <- ggplot() +
  geom_hline(yintercept=0, colour="lightgrey", linetype = 2) +
  geom_vline(xintercept=0, colour="lightgrey", linetype = 2) +
  scale_color_manual(values = col_diet_pred) +
  geom_point(data = hull.uu, aes(x=Axis.1, y=Axis.2, color = diet_pred), alpha = 0.7, size = 4, shape = 16) +
  #stat_ellipse(data = hull.uu,aes(x=Axis.1, y=Axis.2,color=diet_pred),type = "norm")+
  xlab(paste("PCo1 (", round(pcoa.sub.uu$values$Relative_eig[1]*100, 1), "%)")) +
  ylab(paste("PCo2 (", round(pcoa.sub.uu$values$Relative_eig[2]*100, 1), "%)"))  +
  theme_bw() +
  coord_equal() +
  theme(axis.title.x = element_text(size=18, family = "serif"), # remove x-axis labels
        axis.title.y = element_text(size=18, family = "serif"), # remove y-axis labels
        axis.text.x = element_text(size=14, family = "serif"),
        axis.text.y = element_text(size=14, family = "serif"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.title = element_text(family = "serif", size = 16), 
        legend.position = "none",
        #legend.text = element_text(size = 16, family = "serif"),
        #legend.title = element_text(size = 0,family = "serif"),
        strip.text.x = element_text(size=14, family = "serif"))
#+facet_wrap(~ region)
pcoa.diet.pred.uu

pdf(file = "pcoa.diet.pred_elipse.pdf", he = 7, wi = 7)
pcoa.diet.pred.bc
pcoa.diet.pred.wu
pcoa.diet.pred.uu
dev.off()

### PERMANOVA ------
library(vegan)
a.bc = adonis2(otu.bc ~ diet_pred, data = samp_data)
a.wu = adonis2(otu.wu ~ diet_pred, data = samp_data)
a.uu = adonis2(otu.uu ~ diet_pred, data = samp_data)

## *** diet3 *** ----
## __ PCOA plots -----
pcoa.diet3.bc <- ggplot() +
  geom_hline(yintercept=0, colour="lightgrey", linetype = 2) +
  geom_vline(xintercept=0, colour="lightgrey", linetype = 2) +
  scale_color_manual(values = c("darkred","darkgreen" , "darkblue")) +
  geom_point(data = hull.bc, aes(x=Axis.1, y=Axis.2, color = diet3), alpha = 0.7, size = 4, shape = 16) +
  xlab(paste("PCo1 (", round(pcoa.sub.bc$values$Relative_eig[1]*100, 1), "%)")) +
  ylab(paste("PCo2 (", round(pcoa.sub.bc$values$Relative_eig[2]*100, 1), "%)"))  +
  theme_bw() +
  #stat_ellipse(data = hull.bc, aes(x=Axis.1, y=Axis.2,color=diet3),type = "norm")+
  coord_equal() +
  theme(axis.title.x = element_text(size=18, family = "serif"), # remove x-axis labels
        axis.title.y = element_text(size=18, family = "serif"), # remove y-axis labels
        axis.text.x = element_text(size=14, family = "serif"),
        axis.text.y = element_text(size=14, family = "serif"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.title = element_text(family = "serif", size = 16), 
        legend.position = "none",
        #legend.text = element_text(size = 0, family = "serif"),
        #legend.title = element_text(size = 0,family = "serif"),
        strip.text.x = element_text(size=14, family = "serif"))
#+facet_wrap(~ region)
pcoa.diet3.bc

pcoa.diet3.wu <- ggplot() +
  geom_hline(yintercept=0, colour="lightgrey", linetype = 2) +
  geom_vline(xintercept=0, colour="lightgrey", linetype = 2) +
  scale_color_manual(values = c("darkred","darkgreen" , "darkblue")) +
  #stat_ellipse(data = hull.wu, aes(x=Axis.1, y=Axis.2,color=diet3),type = "norm")+
  geom_point(data = hull.wu, aes(x=Axis.1, y=Axis.2, color = diet3), alpha = 0.7, size = 4, shape = 16) +
  xlab(paste("PCo1 (", round(pcoa.sub.wu$values$Relative_eig[1]*100, 1), "%)")) +
  ylab(paste("PCo2 (", round(pcoa.sub.wu$values$Relative_eig[2]*100, 1), "%)"))  +
  theme_bw() +
  coord_equal() +
  theme(axis.title.x = element_text(size=18, family = "serif"), # remove x-axis labels
        axis.title.y = element_text(size=18, family = "serif"), # remove y-axis labels
        axis.text.x = element_text(size=14, family = "serif"),
        axis.text.y = element_text(size=14, family = "serif"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.title = element_text(family = "serif", size = 16), 
        legend.position = "none",
        #legend.text = element_text(size = 0, family = "serif"),
        #legend.title = element_text(size = 0,family = "serif"),
        strip.text.x = element_text(size=14, family = "serif"))
#+facet_wrap(~ region)
pcoa.diet3.wu

pcoa.diet3.uu <- ggplot() +
  geom_hline(yintercept=0, colour="lightgrey", linetype = 2) +
  geom_vline(xintercept=0, colour="lightgrey", linetype = 2) +
  scale_color_manual(values = c("darkred","darkgreen" , "darkblue")) +
  geom_point(data = hull.uu, aes(x=Axis.1, y=Axis.2, color = diet3), alpha = 0.7, size = 4, shape = 16) +
  #stat_ellipse(data = hull.uu,aes(x=Axis.1, y=Axis.2,color=diet3),type = "norm")+
  xlab(paste("PCo1 (", round(pcoa.sub.uu$values$Relative_eig[1]*100, 1), "%)")) +
  ylab(paste("PCo2 (", round(pcoa.sub.uu$values$Relative_eig[2]*100, 1), "%)"))  +
  theme_bw() +
  coord_equal() +
  theme(axis.title.x = element_text(size=18, family = "serif"), # remove x-axis labels
        axis.title.y = element_text(size=18, family = "serif"), # remove y-axis labels
        axis.text.x = element_text(size=14, family = "serif"),
        axis.text.y = element_text(size=14, family = "serif"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.title = element_text(family = "serif", size = 16), 
        legend.position = "none",
        #legend.text = element_text(size = 16, family = "serif"),
        #legend.title = element_text(size = 0,family = "serif"),
        strip.text.x = element_text(size=14, family = "serif"))
#+facet_wrap(~ region)
pcoa.diet3.uu

pdf(file = "pcoa.diet3_elipse.pdf", he = 7, wi = 7)
pcoa.diet3.bc
pcoa.diet3.wu
pcoa.diet3.uu
dev.off()

### PERMANOVA ------
library(vegan)
a.bc.diet3 = adonis(otu.bc ~ family, data = samp_data, strata = samp_data$data2)
a.wu.diet3 = adonis(otu.wu ~ diet3, data = samp_data)
a.uu.diet3 = adonis(otu.uu ~ diet3, data = samp_data)

### CARNIVORES ONLY -----
# Metrics calcul ----
library(phyloseq)
core_carn <- subset_samples(gut_core , diet3 == "Carnivorous")
save(core_carn, file =  "core_carn.RData")
core_rel_carn <- transform_sample_counts(core_carn, function(x) x / sum(x) )
samp_data_carn <- data.frame(sample_data(core_rel_carn))
save(core_rel_carn, file =  "core_rel_carn.RData")

# __ ** Bray-Curtis ----
# PCOA coda & Permanova
library(vegan)
library(ape)
otu.bc.carn <- vegdist(core_rel_carn@otu_table, method = "bray")
pcoa.sub.bc.carn <- pcoa(otu.bc.carn)
pcoa_coord.bc.carn <- pcoa.sub.bc.carn$vectors[,1:3]
# Contruction of the table for graphic 
library(stringr)
hull.bc.carn <- cbind(pcoa_coord.bc.carn, samp_data_carn)
save(hull.bc.carn, file = "hull.bc.carn.RData")
# What is the percentage of the explicative variance? 
library(scales)
paste("Axis 1 :",percent(pcoa.sub.bc.carn$values$Relative_eig[1])) # 12 %
paste("Axis 2 :",percent(pcoa.sub.bc.carn$values$Relative_eig[2])) # 7 %
paste("Axis 3 :",percent(pcoa.sub.bc.carn$values$Relative_eig[3])) # 6 %
paste("Axis 4 :",percent(pcoa.sub.bc.carn$values$Relative_eig[4])) # 5 %
save(otu.bc.carn, hull.bc.carn, file = "beta.bc.carn.RData")

# __ ** Weighted Unifrac ----
# PCOA coda & Permanova
otu.wu.carn <- UniFrac(core_rel_carn, weighted = T, normalized=F, parallel = F, fast=T)
pcoa.sub.wu.carn <- pcoa(otu.wu.carn)
pcoa_coord.wu.carn <- pcoa.sub.wu.carn$vectors[,1:3]
hull.wu.carn <- cbind(pcoa_coord.wu.carn, samp_data_carn)
save(hull.wu.carn, file = "hull.wu.carn.RData")
# What is the percentage of the explicative variance? 
library(scales)
paste("Axis 1 :",percent(pcoa.sub.wu.carn$values$Relative_eig[1])) # 29 %
paste("Axis 2 :",percent(pcoa.sub.wu.carn$values$Relative_eig[2])) # 17 %
paste("Axis 3 :",percent(pcoa.sub.wu.carn$values$Relative_eig[3])) # 10 %
paste("Axis 4 :",percent(pcoa.sub.wu.carn$values$Relative_eig[4])) # 6 %
save(otu.wu.carn, hull.wu.carn, file = "beta.wu.carn.RData")

# __ ** Unweighted Unifrac ----
# PCOA coda & Permanova
otu.uu.carn <- UniFrac(core_rel_carn, weighted = F, normalized=F, parallel = F, fast=T)
pcoa.sub.uu.carn <- pcoa(otu.uu.carn)
pcoa_coord.uu.carn <- pcoa.sub.uu.carn$vectors[,1:3]
hull.uu.carn <- cbind(pcoa_coord.uu.carn, samp_data_carn)
save(hull.uu.carn, file = "hull.uu.carn.RData")

# What is the percentage of the explicative variance? 
library(scales)
paste("Axis 1 :",percent(pcoa.sub.uu.carn$values$Relative_eig[1])) # 14 %
paste("Axis 2 :",percent(pcoa.sub.uu.carn$values$Relative_eig[2])) # 9 %
paste("Axis 3 :",percent(pcoa.sub.uu.carn$values$Relative_eig[3])) # 7 %
paste("Axis 4 :",percent(pcoa.sub.uu.carn$values$Relative_eig[4])) # 5 %
save(otu.uu.carn, hull.uu.carn, file = "beta.uu.carn.RData")


##__ PCOA Plot -----
pcoa.diet.carn.bc <- ggplot() +
  geom_hline(yintercept=0, colour="lightgrey", linetype = 2) +
  geom_vline(xintercept=0, colour="lightgrey", linetype = 2) +
  scale_color_manual(values = col_diet_pred) +
  geom_point(data = hull.bc.carn, aes(x=Axis.1, y=Axis.2, color = diet_pred), alpha = 0.7, size = 4, shape = 16) +
  #stat_ellipse(data = hull.bc.carn, aes(x=Axis.1, y=Axis.2,color=diet_pred),type = "norm")+
  xlab(paste("PCo1 (", round(pcoa.sub.bc.carn$values$Relative_eig[1]*100, 1), "%)")) +
  ylab(paste("PCo2 (", round(pcoa.sub.bc.carn$values$Relative_eig[2]*100, 1), "%)"))  +
  theme_bw() +
  coord_equal() +
  theme(axis.title.x = element_text(size=18, family = "serif"), # remove x-axis labels
        axis.title.y = element_text(size=18, family = "serif"), # remove y-axis labels
        axis.text.x = element_text(size=14, family = "serif"),
        axis.text.y = element_text(size=14, family = "serif"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.title = element_text(family = "serif", size = 16), 
        legend.position = "none",
        #legend.text = element_text(size = 16, family = "serif"),
        #legend.title = element_text(size = 0,family = "serif"),
        strip.text.x = element_text(size=14, family = "serif"))
#+facet_wrap(~ region)
pcoa.diet.carn.bc

pcoa.diet.carn.wu <- ggplot() +
  geom_hline(yintercept=0, colour="lightgrey", linetype = 2) +
  geom_vline(xintercept=0, colour="lightgrey", linetype = 2) +
  scale_color_manual(values = col_diet_pred) +
  geom_point(data = hull.wu.carn, aes(x=Axis.1, y=Axis.2, color = diet_pred), alpha = 0.7, size = 4, shape = 16) +
  #stat_ellipse(data = hull.wu.carn,aes(x=Axis.1, y=Axis.2,color=diet_pred),type = "norm")+
  xlab(paste("PCo1 (", round(pcoa.sub.wu.carn$values$Relative_eig[1]*100, 1), "%)")) +
  ylab(paste("PCo2 (", round(pcoa.sub.wu.carn$values$Relative_eig[2]*100, 1), "%)"))  +
  theme_bw() +
  coord_equal() +
  theme(axis.title.x = element_text(size=18, family = "serif"), # remove x-axis labels
        axis.title.y = element_text(size=18, family = "serif"), # remove y-axis labels
        axis.text.x = element_text(size=14, family = "serif"),
        axis.text.y = element_text(size=14, family = "serif"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.title = element_text(family = "serif", size = 16), 
        legend.position = "none",
        #legend.text = element_text(size = 16, family = "serif"),
        #legend.title = element_text(size = 0,family = "serif"),
        strip.text.x = element_text(size=14, family = "serif"))
#+facet_wrap(~ region)
pcoa.diet.carn.wu

pcoa.diet.carn.uu <- ggplot() +
  geom_hline(yintercept=0, colour="lightgrey", linetype = 2) +
  geom_vline(xintercept=0, colour="lightgrey", linetype = 2) +
  scale_color_manual(values = col_diet_pred) +
  #stat_ellipse(data = hull.uu.carn,aes(x=Axis.1, y=Axis.2,color=diet_pred),type = "norm")+
  geom_point(data = hull.uu.carn, aes(x=Axis.1, y=Axis.2, color = diet_pred), alpha = 0.7, size = 4, shape = 16) +
  xlab(paste("PCo1 (", round(pcoa.sub.uu.carn$values$Relative_eig[1]*100, 1), "%)")) +
  ylab(paste("PCo2 (", round(pcoa.sub.uu.carn$values$Relative_eig[2]*100, 1), "%)"))  +
  theme_bw() +
  coord_equal() +
  theme(axis.title.x = element_text(size=18, family = "serif"), # remove x-axis labels
        axis.title.y = element_text(size=18, family = "serif"), # remove y-axis labels
        axis.text.x = element_text(size=14, family = "serif"),
        axis.text.y = element_text(size=14, family = "serif"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.title = element_text(family = "serif", size = 16), 
        legend.position = "none",
        #legend.text = element_text(size = 16, family = "serif"),
        #legend.title = element_text(size = 0,family = "serif"),
        strip.text.x = element_text(size=14, family = "serif"))
#+ facet_wrap(~ region)
pcoa.diet.carn.uu


pdf(file = "pcoa_carn_elipse.pdf", he = 7, w = 7)
pcoa.diet.carn.bc
pcoa.diet.carn.wu
pcoa.diet.carn.uu
dev.off()

### PERMANOVA ------
library(vegan)
a.bc.carn = adonis2(otu.bc.carn ~ diet_pred, data = samp_data_carn)
a.wu.carn = adonis2(otu.wu.carn ~ diet_pred, data = samp_data_carn)
a.uu.carn = adonis2(otu.uu.carn ~ diet_pred, data = samp_data_carn)

# Save the plot ----
library(gridExtra)
A = pcoa.diet3.bc + xlim(-0.6,0.4) +ylim(-0.5,0.6)
B = pcoa.diet.pred.bc + xlim(-0.6,0.4) +ylim(-0.5,0.6)
C = pcoa.diet.carn.bc + xlim(-0.6,0.4) +ylim(-0.5,0.6)
plot_bc <- grid.arrange(A, B, C, ncol=3, nrow = 1)

D = pcoa.diet3.wu + xlim(-0.5,0.6) +ylim(-0.4,0.8)
E = pcoa.diet.pred.wu + xlim(-0.5,0.6) +ylim(-0.4,0.8)
F = pcoa.diet.carn.wu + xlim(-0.5,0.6) +ylim(-0.4,0.8)
plot_wu <- grid.arrange(D, E, F, ncol=3, nrow = 1)

G = pcoa.diet3.uu + xlim(-0.4,0.4) +ylim(-0.6,0.3)
H = pcoa.diet.pred.uu + xlim(-0.4,0.4) +ylim(-0.6,0.3)
I = pcoa.diet.carn.uu + xlim(-0.4,0.4) +ylim(-0.6,0.3)
plot_uu <- grid.arrange(G, H, I, ncol=3, nrow = 1)

pdf(file = "pcoa.plots.pdf", he = 7, wi = 10)
plot_bc
plot_wu
plot_uu
dev.off()



### Venn diagram ----
source("http://raw.github.com/nielshanson/mp_tutorial/master/downstream_analysis_r/code/venn_diagram3.r")
library(phyloseq)

merged = merge_samples(gut_core, "diet3") #merge samples for herbivores, carnivores and the omnivores
rownames(merged@otu_table@.Data)
min(rowSums(merged@otu_table@.Data)) #how many reads per sample
set.seed(10000)
venn = rarefy_even_depth(merged, sample.size = min(rowSums(merged@otu_table@.Data)))

herb <- colnames(venn@otu_table@.Data)[venn@otu_table@.Data["Herbivorous",] > 0]
carn <- colnames(venn@otu_table@.Data)[venn@otu_table@.Data["Carnivorous",] > 0]
omni <- colnames(venn@otu_table@.Data)[venn@otu_table@.Data["Omnivorous",] > 0]

pdf(file = "venn_diet.pdf", he = 7, wi = 7)
venn_diagram3(herb,carn,omni, "Herbivores", "Carnivores", "Omnivorous", colors= c("darkgreen","darkred","darkblue"), euler=FALSE)
dev.off()


## Betadispersion ----
samp_data_diet <- core_rel %>% sample_data() %>% as.data.frame()
disp_diet.bc <- betadisper(otu.bc, samp_data_diet$diet_pred)
Distances.bc <- disp_diet.bc$distances
disp_diet.wu <- betadisper(otu.wu, samp_data_diet$diet_pred)
Distances.wu <- disp_diet.wu$distances
sig.diet.wu <- permutest(disp_diet.wu)
disp_diet.uu <- betadisper(otu.uu, samp_data_diet$diet_pred)
Distances.uu <- disp_diet.uu$distances
sig.diet.uu <- permutest(disp_diet.uu)
samp_data_diet <- cbind("ID"= rownames(samp_data_diet), samp_data_diet)
disp.diet <- samp_data_diet %>% select(ID, tax1, diet_pred)
Distances <- rbind(cbind(disp.diet , "Betadisp" = Distances.bc, "Index" = rep("Bray-Curtis")), 
                   cbind(disp.diet , "Betadisp" = Distances.wu, "Index" = rep("W-Uni")),
                   cbind(disp.diet, "Betadisp" = Distances.uu, "Index" = rep("U-Uni")))

Distances$diet_pred <- factor(Distances$diet_pred , levels = c("Herbivores Microvores Detritivores",
                                                               "Planktivores",
                                                               "Corallivores",
                                                               "sessile invertivores",
                                                               "Microinvertivores",
                                                               "Crustacivores",
                                                               "Macroinvertivores",
                                                               "Piscivores"
                                                               ))

write.table(Distances , file ="Distances_betadisp_diet.txt", sep ="\t", row.names = F)
col_diet_pred <- c('Corallivores' =  'indianred1',
                   'Crustacivores' = 'orange',
                   'Herbivores Microvores Detritivores' = 'darkgreen',
                   'Macroinvertivores' = 'darkslategray',
                   'Microinvertivores' = 'bisque3',
                   'Piscivores' = 'darkblue',
                   'Planktivores' = 'darkred',
                   'sessile invertivores' = 'coral4',
                   'Omnivorous' = 'gray',
                   'Carnivorous' = 'darkslateblue')

plot.disp.diet =  ggplot(Distances, aes(x = Betadisp , y = diet_pred, color = diet_pred)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  geom_boxplot(alpha=0.1, outlier.colour = NA) +
  scale_color_manual(values = col_diet_pred)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=0, family = "serif", face = "italic"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=12, family = "serif"),
        legend.position="none",
        legend.text = element_text(family = "serif"),
        legend.title = element_text(family = "serif")) +
  facet_wrap( ~ Index, nrow=1, ncol=3, scales = "free")

pdf(file ="plot.disp.diet.pdf", he = 7 , wi = 7)
plot.disp.diet
dev.off()

library(dunn.test)
kruskal.test(Distances$Betadisp[Distances$Index == "W-Uni"] , Distances$diet_pred[Distances$Index == "W-Uni"])
dunn.test(Distances$Betadisp[Distances$Index == "Bray-Curtis"] , Distances$diet_pred[Distances$Index == "Bray-Curtis"], method = "bonferroni")

### Correlation ----
trophic <- read.table(file = "../Alpha/box_alpha_troph.txt", sep = "\t", header = T)
trophic <- trophic %>% select(tax1, trophic_position) %>% unique()

Distances <- left_join(Distances, trophic, by = "tax1")

bc_cor <- cor.test(Distances$Betadisp[Distances$Index == "Bray-Curtis"], Distances$trophic_position[Distances$Index == "Bray-Curtis"], method=c("spearman"))
wu_cor <- cor.test(Distances$Betadisp[Distances$Index == "W-Uni"], Distances$trophic_position[Distances$Index == "W-Uni"], method=c("spearman"))
uu_cor <- cor.test(Distances$Betadisp[Distances$Index == "U-Uni"], Distances$trophic_position[Distances$Index == "U-Uni"], method=c("spearman"))

correlation_beta_diet <- ggplot(Distances, aes(x=trophic_position , y= Betadisp)) + 
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm")+
  xlab("Trophic position")+
  ylab("Beta dispersion")+
  theme(axis.title.x = element_text(family = "serif",size = 16),
        axis.title.y = element_text(family = "serif",size = 16),
        axis.title = element_text(family = "serif",size = 14),
        axis.text.y = element_text(size=14, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 14, angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=14, family = "serif"))+
  facet_wrap(~ Index, scales  = "free")

pdf(file = 'correlation_diet_beta.pdf', he = 5, wi = 7)
correlation_beta_diet
dev.off()


## Unique genera ----
# Occurrence
length(which(sample_sums(subset_taxa(gut_core, Genus == "NS10_marine_group")) > 0))
# Mean abudance in reservoirs 
mean(sample_sums(subset_taxa(gut_core, Genus == "NS10_marine_group"))[names(sample_sums(subset_taxa(gut_core, Genus == "NS10_marine_group"))) %in% 
        names(which(sample_sums(subset_taxa(gut_core, Genus == "NS10_marine_group")) > 0))]/
  sample_sums(gut_core)[names(sample_sums(subset_taxa(gut_core, Genus == "NS10_marine_group"))) %in% names(which(sample_sums(subset_taxa(gut_core, Genus == "NS10_marine_group")) > 0))] * 100)

samp_data[rownames(samp_data) %in% names(which(sample_sums(subset_taxa(gut_core, Genus == "NS10_marine_group")) > 0)),]

sd(sample_sums(subset_taxa(gut_core, Genus == "NS10_marine_group"))[names(sample_sums(subset_taxa(gut_core, Genus == "NS10_marine_group"))) %in% 
                                                              names(which(sample_sums(subset_taxa(gut_core, Genus == "NS10_marine_group")) > 0))]/
       sample_sums(gut_core)[names(sample_sums(subset_taxa(gut_core, Genus == "NS10_marine_group"))) %in% names(which(sample_sums(subset_taxa(gut_core, Genus == "NS10_marine_group")) > 0))] * 100)

chitin.tax <- subset_taxa(rff_core , Genus %in% c("Vibrio" , "Photobacterium"))
sum(sample_sums(subset_samples(chitin.tax, diet3 == "Herbivorous")))/sum(sample_sums(subset_samples(rff_core, diet3 == "Herbivorous"))) *100
sum(sample_sums(subset_samples(chitin.tax, diet3 == "Carnivorous")))/sum(sample_sums(subset_samples(rff_core, diet3 == "Carnivorous"))) *100
