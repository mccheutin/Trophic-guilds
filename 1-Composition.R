## 00 - Setup assignation and Composition ------

## Trophic guilds assignation -----
setwd("/Users/marie-charlottecheutin/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Carnivory")
trophic_guild <- read.table(file = "trophic_guilds.txt", header = T, sep = "\t")
head(trophic_guild)
colnames(trophic_guild)[2] <- "tax2" 

env <- as.data.frame(sample_data(gut_core))
write.table(env, file = "physeq_object/env.txt", row.names = F, sep ="\t", quote = F)
env <- read.table(file = "physeq_object/env.txt", header = T, sep = "\t")

env_tg <-  left_join(env , trophic_guild, by = "tax2")
write.table(env_tg, file = "physeq_object/env_tg.txt", row.names = F, sep ="\t", quote = F)
env_tg <- read.table(file = "physeq_object/env_tg.txt", header = T, sep = "\t")
rownames(env_tg) <- sample_names(gut_core)
sample_data(gut_core) <- env_tg
save(gut_core, file = "physeq_object/gut_core.RData")


# __ Figure Description ----
core_physeq_genus <- gut_core %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at order level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Genus)

core_physeq_genus$Phylum = as.character(core_physeq_genus$Phylum) # Avoid error message with factor for next step
sort(table(factor(core_physeq_genus$Phylum)), T)
core_physeq_genus[!core_physeq_genus$Phylum %in% 
                    c("Proteobacteria","Bacteroidetes","Cyanobacteria","Firmicutes","Actinobacteria",
                      "Planctomycetes", "Spirochaetes", "Verrucomicrobia","Fusobacteria","Tenericutes"),
                  which(names(core_physeq_genus) == "Phylum", T)] <- "Other" #Change phylum to other for those not included in the list

core_physeq_genus$Class = as.character(core_physeq_genus$Class) # Avoid error message with factor for next step
sort(table(factor(core_physeq_genus$Class)), T)
core_physeq_genus[!core_physeq_genus$Class %in% 
                    rownames(which(sort(table(factor(core_physeq_genus$Class)), T) > 1000,T)),
                  which(names(core_physeq_genus) == "Class", T)] <- "Other" 

save(core_physeq_genus, file = "physeq_object/core_physeq_genus.RData")

# __ Figure Barplot  ----
core_physeq_genus$diet_pred <- factor(core_physeq_genus$diet_pred , levels = c("Piscivores",
                                                                               "Macroinvertivores", 
                                                                               "Crustacivores",
                                                                               "Microinvertivores",
                                                                               "sessile invertivores",
                                                                               "Corallivores",
                                                                               "Planktivores",
                                                                               "Herbivores Microvores Detritivores"))

levels(core_physeq_genus$diet_pred)[5] <- "Sessile invertivores"
core_physeq_genus$Phylum <- factor(core_physeq_genus$Phylum , levels = c(names(sort(table(factor(core_physeq_genus$Phylum)), T))[names(sort(table(factor(core_physeq_genus$Phylum)), T)) != "Other"],
                                                                         "Other"))


save(core_physeq_genus , file ="physeq_object/core_physeq_genus_modif.RData")

## Phylum plot ----
phylum_colors <- c('Proteobacteria'='aquamarine4',
                   'Firmicutes'='brown4',
                   'Bacteroidetes'='cornflowerblue',
                   'Cyanobacteria'='#DE6554',
                   'Actinobacteria'='darkblue',
                   'Planctomycetes'='darkslategray',
                   'Verrucomicrobia'='darkgreen',
                   'Tenericutes'='#CA8EA7',
                   'Spirochaetes'='darkolivegreen3',
                   'Fusobacteria'='orange',
                   'Other'= 'black')

sum_phylum <- core_physeq_genus %>% group_by(Phylum,diet_pred) %>% summarise(Sum = sum(Abundance))

plot_phylum <- ggplot(sum_phylum, aes(x = Sum, y = diet_pred, fill = Phylum)) + 
  geom_bar(stat="identity", position="fill") + 
  xlab("Relative Abundance") +
  scale_x_continuous(expand = c(0,0)) + #remove the space below the 0 of the y axis in the graph
  ggtitle("Phylum") +
  scale_fill_manual(values= phylum_colors) +
  theme_bw() +
  theme(axis.line = element_line(size = 0), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(), 
        axis.title.x = element_text(family = "serif",size = 15),
        axis.title.y = element_text(family = "serif",size = 0),
        axis.text = element_text(size = 14), 
        axis.text.y = element_text(size = 15, family = "serif"), 
        axis.text.x = element_text(family = "serif", size = 14), 
        plot.title = element_text(family = "serif", size = 20), 
        legend.text = element_text(size = 15,  family = "serif"),
        legend.title = element_text(size = 0, family = "serif"))

plot_phylum

### Class plot -----
core_physeq_genus$Class <- factor(core_physeq_genus$Class , levels = c(names(sort(table(factor(core_physeq_genus$Class)), T))[names(sort(table(factor(core_physeq_genus$Class)), T)) != "Other"],
                                                                       "Other"))

class_colors <- c('Alphaproteobacteria'='lightseagreen',
                  'Clostridia'='orange3',
                  'Gammaproteobacteria'='aquamarine4',
                  'Bacteroidia'='cornflowerblue',
                  'Oxyphotobacteria'='#DE6554',
                  'Deltaproteobacteria'='skyblue4',
                  'Bacilli'='brown4',
                  'Planctomycetacia'='darkslategray',
                  'Verrucomicrobiae'='darkgreen',
                  'Actinobacteria'='darkblue',
                  'Erysipelotrichia'='yellow',
                  'Mollicutes'='#CA8EA7',
                  'Fusobacteriia'='orange',
                  'Spirochaetia'='darkolivegreen3',
                  'Other' = 'black')

sum_class <- core_physeq_genus %>% group_by(Class,diet_pred) %>% summarise(Sum = sum(Abundance))

plot_class <- ggplot(sum_class, aes(x = Sum, y = diet_pred, fill = Class)) + 
  geom_bar(stat="identity", position="fill") + 
  xlab("Relative Abundance") +
  scale_x_continuous(expand = c(0,0)) + #remove the space below the 0 of the y axis in the graph
  ggtitle("Class") +
  scale_fill_manual(values= class_colors) +
  theme_bw() +
  theme(axis.line = element_line(size = 0), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(), 
        axis.title.x = element_text(family = "serif",size = 15), 
        axis.title.y = element_text(family = "serif",size = 0), 
        axis.text = element_text(size = 14), 
        axis.text.y = element_text(size = 15, family = "serif"), 
        axis.text.x = element_text(family = "serif", size = 13), 
        plot.title = element_text(family = "serif", size = 20), 
        legend.text = element_text(size = 15,  family = "serif"),
        legend.title = element_text(size = 0, family = "serif"))

plot_class

## Genus plot ----
load("physeq_object/core_physeq_genus_modif.RData")
Genus_Ab <- core_physeq_genus %>% 
  group_by(Genus) %>%
  summarize(Sum = sum(Abundance, na.rm=TRUE))


Genus_Ab <- Genus_Ab[order(-Genus_Ab$Sum),]
Genus_Ab$Sum <- Genus_Ab$Sum/sum(Genus_Ab$Sum)
Genus_15 <- as.character(Genus_Ab[1:15,]$Genus)
Genus_other <- Genus_Ab %>% filter(!Genus %in% Genus_15) %>% select(Genus)
Genus_other <- as.character(Genus_other$Genus)

core_physeq_genus$Genus = as.character(core_physeq_genus$Genus) # Avoid error message with factor for next step
core_physeq_genus[!core_physeq_genus$Genus %in% Genus_15,
                  which(names(core_physeq_genus) == "Genus", T)] <- "Other" #Change phylum to other for those not included in the list

core_physeq_genus$Genus <- factor(core_physeq_genus$Genus , levels = c(Genus_15, "Other"))
n <- length(levels(as.factor(core_physeq_genus$Genus)))
genus_colors <- c(distinctColorPalette(n-1),"black")

sum_genus <- core_physeq_genus %>% group_by(Genus,diet_pred) %>% summarise(Sum = sum(Abundance))
plot_genus <- ggplot(sum_genus, aes(x = Sum, y = diet_pred, fill = Genus)) +
  geom_bar(stat="identity", position="fill") + 
  xlab("Relative Abundance") +
  scale_x_continuous(expand = c(0,0)) + #remove the space below the 0 of the y axis in the graph
  ggtitle("Genus") +
  scale_fill_manual(values= genus_colors) +
  theme_bw() +
  theme(axis.line = element_line(size = 0), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(), 
        axis.title.x = element_text(family = "serif",size = 15), 
        axis.title.y = element_text(family = "serif",size = 0), 
        axis.text = element_text(size = 14), 
        axis.text.y = element_text(size = 15, family = "serif"), 
        axis.text.x = element_text(family = "serif", size = 13), 
        plot.title = element_text(family = "serif", size = 20), 
        legend.text = element_text(size = 15,  family = "serif"),
        legend.title = element_text(size = 0, family = "serif"))

plot_genus

pdf(file ="compo_files.pdf" , he = 5 , wi = 10)
plot_phylum
plot_class
plot_genus
dev.off()

##### Merging diet_pred ----
merged_core = merge_samples(gut_core, "diet_pred")
diet_physeq_genus <- merged_core %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at order level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Genus)


cora_gen <- diet_physeq_genus %>% filter(Sample == "Corallivores") %>% group_by(Genus)
cora_gen <- cora_gen[cora_gen$Abundance >0, ]$Genus
macroinv_gen <- diet_physeq_genus %>% filter(Sample == "Macroinvertivores") %>% group_by(Genus)
macroinv_gen <- macroinv_gen[macroinv_gen$Abundance >0, ]$Genus
plank_gen <- diet_physeq_genus %>% filter(Sample == "Planktivores") %>% group_by(Genus)
plank_gen <- plank_gen[plank_gen$Abundance >0, ]$Genus
crust_gen <- diet_physeq_genus %>% filter(Sample == "Crustacivores") %>% group_by(Genus)
crust_gen <- crust_gen[crust_gen$Abundance >0, ]$Genus
microinv_gen <- diet_physeq_genus %>% filter(Sample == "Microinvertivores") %>% group_by(Genus)
microinv_gen <- microinv_gen[microinv_gen$Abundance >0, ]$Genus
si_gen <- diet_physeq_genus %>% filter(Sample == "sessile invertivores") %>% group_by(Genus)
si_gen <- si_gen[si_gen$Abundance >0, ]$Genus
pisci_gen <- diet_physeq_genus %>% filter(Sample == "Piscivores") %>% group_by(Genus)
pisci_gen <- pisci_gen[pisci_gen$Abundance >0, ]$Genus
hmd_gen <- diet_physeq_genus %>% filter(Sample == "Herbivores Microvores Detritivores") %>% group_by(Genus)
hmd_gen <- hmd_gen[hmd_gen$Abundance >0, ]$Genus

table_genus <- as.data.frame(rbind(cbind("Genus" = cora_gen ,"Diet" = rep("Corallivores")),
                                   cbind("Genus" = macroinv_gen ,"Diet" = rep("Macroinvertivores")),
                                   cbind("Genus" = plank_gen ,"Diet" = rep("Planktivores")),
                                   cbind("Genus" = crust_gen ,"Diet" = rep("Crustacivores")),
                                   cbind("Genus" = microinv_gen ,"Diet" = rep("Microinvertivores")),
                                   cbind("Genus" = si_gen ,"Diet" = rep("SI")),
                                   cbind("Genus" = pisci_gen ,"Diet" = rep("Piscivores")),
                                   cbind("Genus" = hmd_gen ,"Diet" = rep("HMD"))))

unique_genes <- names(which(table(table_genus$Genus) == 1))
unique_table <- table_genus %>% filter(Genus %in% unique_genes)
write.table(unique_table, file = "unique_genera_per_diet.txt", sep ="\t", row.names = F, quote = F)

diet_physeq_genus %>% filter(Genus %in% unique_genes) %>% summarize(Sum = sum(Abundance, na.rm=TRUE))

## Lokking for the abundances of Genera in class of diets ----
diet3_core = merge_samples(gut_core, "diet3")

diet3_physeq_genus <- diet3_core %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at order level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Genus)

herb_physeq <- diet3_physeq_genus %>% filter(Sample == "Herbivorous")  
herb_physeq %>% filter(Order== "Nostocales") %>% select(Abundance) %>% sum()
carn_physeq <- diet3_physeq_genus %>% filter(Sample == "Carnivorous") 
omni_physeq <- diet3_physeq_genus %>% filter(Sample == "Omnivorous") 


