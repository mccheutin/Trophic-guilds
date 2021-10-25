### How diversity is influenced by the diet?
### 
### Alpha diversity -----
setwd("/Users/marie-charlottecheutin/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Carnivory")
dir.create("./Alpha", recursive = T)
setwd("./Alpha")

load("../physeq_object/gut_core.RData")
sort(sample_sums(gut_core))

rff_core <- rarefy_even_depth(gut_core, sample.size = 1060)
rff_core

goods <-
  function(com){
    no.seqs <- rowSums(com)
    sing <- com==1
    no.sing <- apply(sing, 1, sum)
    goods <- 100*(1-no.sing/no.seqs)
    goods.sum <- cbind(no.sing, no.seqs, goods)
    goods.sum <- as.data.frame(goods.sum)
    return(goods.sum)
  }

cov_core_rff <- goods(rff_core@otu_table@.Data)
paste0(mean(cov_core_rff$goods)," +- " ,sd(cov_core_rff$goods))
save(rff_core, file = "rff_core.RData")

# __ Metrics calcul ----
library(picante)
library(entropart)
box1 = estimate_richness(rff_core , measures = c("Shannon", "Chao1"))
box1$Shannon = exp(box1$Shannon)
box2 = pd(samp= rff_core@otu_table , tree = rff_core@phy_tree , include.root =  T)
box3 <- as.data.frame(cbind("Sample" = rownames(box1), box1, box2))
db <- rff_core %>% sample_data() %>% as.data.frame()
db <- cbind("Sample" = rownames(db), db)
library(stringr)
db$Sample <- str_replace_all(db$Sample, "-", ".")
box3 <- left_join(box3, db , by = "Sample") 
save(box3, file = "box.alpha.RData")

# *** A) Alpha diet ----
box3_obs <- box3$SR
box3_PD <- box3$PD
box3_shannon <- box3$Shannon
box3_chao <- box3$Chao1
box3_value <- as.data.frame(rbind(cbind("variable" = rep("Observed"), "value" = box3_obs),
                                      cbind("variable" = rep("PD"), "value" = box3_PD),
                                      cbind("variable" = rep("Shannon"), "value" = box3_shannon),
                                      cbind("variable" = rep("Chao1"), "value" = box3_chao)))

box_plot <- box3 %>% select(Sample, tax1, tax2, family, diet3,diet_pred)
box3_plot <- cbind(box_plot, box3_value)
box3_plot$value <- as.numeric(box3_plot$value)
box3_plot$variable <- factor(box3_plot$variable , levels = c("Observed", "PD","Shannon","Chao1"))
save(box3_plot, file = "box3_plot.RData")
load("box3_plot.RData")

box3_plot$diet_pred <- str_replace_all(box3_plot$diet_pred, c("Herbivores Microvores Detritivores" = "HMD",
                                                   "sessile invertivores" = "Sessile invertivores"))
box3_plot$diet_pred <- factor(box3_plot$diet_pred , levels = rev(c("Piscivores",
                                                               "Macroinvertivores",
                                                               "Crustacivores",
                                                               "Microinvertivores",
                                                               "Sessile invertivores",
                                                               "Corallivores",
                                                               "Planktivores",
                                                               "HMD")))

col_diet_pred <- c('Corallivores' =  'indianred1',
                   'Crustacivores' = 'orange',
                   'HMD' = 'darkgreen',
                   'Macroinvertivores' = 'darkslategray',
                   'Microinvertivores' = 'bisque3',
                   'Piscivores' = 'darkblue',
                   'Planktivores' = 'darkred',
                   'Sessile invertivores' = 'coral4')



diet_alpha_div =  ggplot(box3_plot %>% filter(variable %in% c("Shannon", "Chao1")), aes(x = as.numeric(value) , y = diet_pred, color = diet_pred)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_colour_manual(values= col_diet_pred) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=14, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 14, angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=14, family = "serif"),
        legend.position="none")+
  facet_grid( ~ variable,  scales = "free")
diet_alpha_div

pdf(file = "alpha_plot_diet.V2.pdf" , he = 5 , wi = 7)
diet_alpha_div
dev.off()

# __ * Dunn test ----
box3_plot <- load("box3_plot.RData")
library(dunn.test)
obs.diet.pred <- dunn.test(as.numeric(box3_plot$value)[box3_plot$variable == "Observed"] , box3_plot$diet_pred[box3_plot$variable=="Observed"], method = "bonferroni")
shannon.diet.pred <- dunn.test(as.numeric(box3_plot$value)[box3_plot$variable == "Shannon"] , box3_plot$diet_pred[box3_plot$variable=="Shannon"], method = "bonferroni")
pd.diet.pred <- dunn.test(as.numeric(box3_plot$value)[box3_plot$variable == "PD"] , box3_plot$diet_pred[box3_plot$variable=="PD"], method = "bonferroni")
chao.diet.pred <- dunn.test(as.numeric(box3_plot$value)[box3_plot$variable == "Chao1"] , box3_plot$diet_pred[box3_plot$variable=="Chao1"], method = "bonferroni")

comp.diet.pred <- as.data.frame(cbind(obs.diet.pred$comparisons, 
                                      round(obs.diet.pred$P.adjusted,2) ,
                                      round(shannon.diet.pred$P.adjusted,2), 
                                      round(pd.diet.pred$P.adjusted,2),
                                      round(chao.diet.pred$P.adjusted,2)))

colnames(comp.diet.pred) <- c("Comparisons",  "Obs", "Shannon", "PD", "Chao1")
write.table(comp.diet.pred, file = "comp.diet.pred.txt", sep = "\t", quote = F, row.names = F)

obs.diet3 <- dunn.test(as.numeric(box3_plot$value)[box3_plot$variable == "Observed"] , box3_plot$diet3[box3_plot$variable=="Observed"], method = "bonferroni")
shannon.diet3 <- dunn.test(as.numeric(box3_plot$value)[box3_plot$variable == "Shannon"] , box3_plot$diet3[box3_plot$variable=="Shannon"], method = "bonferroni")
pd.diet3 <- dunn.test(as.numeric(box3_plot$value)[box3_plot$variable == "PD"] , box3_plot$diet3[box3_plot$variable=="PD"], method = "bonferroni")
chao.diet3 <- dunn.test(as.numeric(box3_plot$value)[box3_plot$variable == "Chao1"] , box3_plot$diet3[box3_plot$variable=="Chao1"], method = "bonferroni")


comp.diet3 <- as.data.frame(cbind(obs.diet3$comparisons, 
                                  round(obs.diet3$P.adjusted,2) , 
                                  round(shannon.diet3$P.adjusted,2), 
                                  round(pd.diet3$P.adjusted,2),
                                  round(chao.diet3$P.adjusted,2)))

colnames(comp.diet3) <- c("Comparisons",  "Obs", "Shannon", "PD", "Chao1")
write.table(comp.diet3, file = "comp.diet3.txt", sep = "\t", quote = F, row.names = F)



### Globally, it seems that richness of herbivores > omnivorous >= carnivorous
box3_plot <- read.table("box3_plot.txt", header = T, sep = "\t")
box3_plot$value <- as.numeric(box3_plot$value)

library(dplyr)
mean_sd_diet_pred <- box3_plot %>% 
  select(c(diet3,diet_pred, variable, value)) %>%
  filter(variable == "Chao1") %>%
  group_by(diet_pred) %>%
  summarize(Mean = mean(value, na.rm=TRUE), SD = sd(value, na.rm = T))

mean_sd_diet3 <- box3_plot %>% 
  select(c(diet3,diet_pred, variable, value)) %>%
  filter(variable == "Observed") %>%
  group_by(diet3) %>%
  summarize(Mean = mean(value, na.rm=TRUE), SD = sd(value, na.rm = T))

mean_sd_position <- box3_plot %>% 
  select(c(diet3,diet_pred,trophic_position, variable, value)) %>%
  filter(variable == "Richness observed") %>%
  group_by(trophic_position) %>%
  summarize(Mean = mean(value, na.rm=TRUE), SD = sd(value, na.rm = T))

box3_plot %>% 
  select(c(diet3,diet_pred,trophic_position, variable, value)) %>%
  filter(variable == "Richness observed")

#### Correlation between richness and trophic position ----
troph_dat <- read.table("troph.txt", sep ="\t", header = T)
troph_dat
box_alpha_troph <- left_join(box3_plot, troph_dat[,c(1,7)],by = "tax1")
box_alpha_troph <- unique(box_alpha_troph)
write.table(box_alpha_troph, file = "box_alpha_troph.txt" , row.names = F, sep = "\t")


library("ggpubr")
ggqqplot(box_alpha_troph$value[box_alpha_troph$variable == "Observed"])
shapiro.test(box_alpha_troph[box_alpha_troph$variable == "Chao1",]$value)

library(ggplot2)
box_alpha_troph$variable <- factor(box_alpha_troph$variable, levels = c("Observed", "Shannon", "PD", "Chao1"))
correlation_alpha_diet <- ggplot(box_alpha_troph %>% filter(variable %in% c("Shannon", "Chao1")), aes(x= as.numeric(trophic_position) , y= as.numeric(value))) + 
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm")+
  xlab("Trophic position")+
  ylab("Value")+
  theme(axis.title.x = element_text(family = "serif",size = 16),
        axis.title.y = element_text(family = "serif",size = 16),
        axis.title = element_text(family = "serif",size = 14),
        axis.text.y = element_text(size=14, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 14, angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=14, family = "serif"))+
  facet_wrap(~ variable, scales  = "free")

pdf(file = 'correlation_diet_sha_cha.pdf', he = 5, wi = 7)
correlation_alpha_diet
dev.off()

richness <- cor.test(as.numeric(box_alpha_troph[box_alpha_troph$variable == "Observed",]$trophic_position), as.numeric(box_alpha_troph[box_alpha_troph$variable == "Observed",]$value), method=c("spearman"))
shannon <- cor.test(as.numeric(box_alpha_troph[box_alpha_troph$variable == "Shannon",]$trophic_position), as.numeric(box_alpha_troph[box_alpha_troph$variable == "Shannon",]$value), method=c("spearman"))
pd <- cor.test(as.numeric(box_alpha_troph[box_alpha_troph$variable == "PD",]$trophic_position), as.numeric(box_alpha_troph[box_alpha_troph$variable == "PD",]$value), method=c("spearman"))
chao <- cor.test(as.numeric(box_alpha_troph[box_alpha_troph$variable == "Chao1",]$trophic_position), as.numeric(box_alpha_troph[box_alpha_troph$variable == "Chao1",]$value), method=c("spearman"))

cor.val.richness <- round(cor(as.numeric(box_alpha_troph[box_alpha_troph$variable == "Observed",]$value), 
                              as.numeric(box_alpha_troph[box_alpha_troph$variable == "Observed",]$trophic_position), 
                              method = "spearman"),
                          2)

cor.val.shannon <- round(cor(as.numeric(box_alpha_troph[box_alpha_troph$variable == "Shannon",]$value), 
                             as.numeric(box_alpha_troph[box_alpha_troph$variable == "Shannon",]$trophic_position), 
                             method = "spearman"),
                         2)

cor.val.pd <- round(cor(as.numeric(box_alpha_troph[box_alpha_troph$variable == "PD",]$value), 
                        as.numeric(box_alpha_troph[box_alpha_troph$variable == "PD",]$trophic_position), 
                        method = "spearman"),
                    2)

cor.val.chao <- round(cor(as.numeric(box_alpha_troph[box_alpha_troph$variable == "Chao1",]$value), 
                        as.numeric(box_alpha_troph[box_alpha_troph$variable == "Chao1",]$trophic_position), 
                        method = "spearman"),
                    2)
