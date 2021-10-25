### Core Descritpion ------------------------------------------
#Setup ----
knitr::opts_chunk$set(eval = FALSE)
remove(list = ls())
cran_packages   <- c("knitr", "phyloseqGraphTest", "phyloseq", "shiny", "microbiome",
                     "tidyverse", "miniUI", "caret", "pls", "e1071", "ggplot2", 
                     "randomForest","entropart", "vegan", "plyr", "dplyr", "here",
                     "ggrepel", "nlme", "R.utils", "gridExtra","grid", "googledrive", 
                     "googlesheets", "phangorn", "devtools", "rmarkdown", "sys",
                     "reshape2", "devtools", "PMA","structSSI","ade4", "ape",
                     "Biostrings", "igraph", "ggnetwork", "intergraph", "ips",
                     "scales", "kableExtra", "pgirmess", "treemap", "knitr","kableExtra",
                     "rstudioapi" ,"data.table","DT","pander","formatR","grDevices","svgPanZoom",
                     "RCurl","plotly","pairwiseAdonis", "stringr")
github_packages <- c("jfukuyama/phyloseqGraphTest")
bioc_packages   <- c("phyloseq", "genefilter", "impute", "dada2", "DECIPHER")
# Install CRAN packages (if not already installed)
#Some packages would be not availbale for your R version
inst <- cran_packages %in% installed.packages()
if (any(! inst)) {
  install.packages(cran_packages[!inst], repos = "http://cran.rstudio.com/") }
# 
inst <- github_packages %in% installed.packages()
if (any(! inst)) {
  devtools::install_github(github_packages[!inst]) }

# Load libraries
sapply(c(cran_packages, bioc_packages), require, character.only = TRUE)
sessionInfo()
set.seed(1000)

setwd("~/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Carnivory")

#Data preparation ----
load("/Users/marie-charlottecheutin/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Meta_analisis/Physeq_objects/gut_core.RData")

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

save(core_physeq_genus, file = "physeq_object/core_physeq_genus.RData")

core_physeq_genus$Class = as.character(core_physeq_genus$Class) # Avoid error message with factor for next step
sort(table(factor(core_physeq_genus$Class)), T)
core_physeq_genus[!core_physeq_genus$Class %in% 
                    rownames(which(sort(table(factor(core_physeq_genus$Class)), T) > 1000,T)),
                  which(names(core_physeq_genus) == "Class", T)] <- "Other" 

core_physeq_genus_rel_herb <- core_physeq_genus %>%
  select(c(Abundance, diet3, diet_pred, Phylum, Class, Order, Family, Genus)) %>%
  filter(diet3 == "Herbivorous")
core_physeq_genus_rel_herb$Abundance <- core_physeq_genus_rel_herb$Abundance/sum(core_physeq_genus_rel_herb$Abundance)

core_physeq_genus_rel_carn <- core_physeq_genus %>%
  select(c(Abundance, diet3, diet_pred, Phylum, Class, Order, Family, Genus)) %>%
  filter(diet3 == "Carnivorous")
core_physeq_genus_rel_carn$Abundance <- core_physeq_genus_rel_carn$Abundance/sum(core_physeq_genus_rel_carn$Abundance)

core_physeq_genus_rel_omni <- core_physeq_genus %>%
  select(c(Abundance, diet3, diet_pred, Phylum, Class, Order, Family, Genus)) %>%
  filter(diet3 == "Omnivorous")
core_physeq_genus_rel_omni$Abundance <- core_physeq_genus_rel_omni$Abundance/sum(core_physeq_genus_rel_omni$Abundance)

save(core_physeq_genus_rel_carn, core_physeq_genus_rel_herb, core_physeq_genus_rel_omni, file = "physeq_object/core_physeq_rel.RData")

core_physeq_genus_rel_herb %>%
  group_by(Class) %>%
  summarize(Sum = sum(Abundance, na.rm=TRUE)*100)

core_physeq_genus_rel_carn %>%
  group_by(Class) %>%
  summarize(Sum = sum(Abundance, na.rm=TRUE)*100)

core_physeq_genus_rel_omni %>%
  group_by(Class) %>%
  summarize(Sum = sum(Abundance, na.rm=TRUE)*100)

# __ Figure Barplot  ----
phylum_colors <- c('Acidobacteria'='lavenderblush4',
                   'Actinobacteria'='darkblue',
                   'Armatimonadetes'='cadetblue3',
                   'Bacteroidetes'='cornflowerblue',
                   'Calditrichaeota'='azure3',
                   'Chloroflexi'='#DCE1D2',
                   'Cyanobacteria'='#DE6554',
                   'Dadabacteria'='brown2',
                   'Deferribacteres'='darkslategray1',
                   'Deinococcus-Thermus'='salmon4',
                   'Dependentiae'='sandybrown',
                   'Epsilonbacteraeota'='darkslateblue',
                   'Elusimicrobia'='plum1',
                   'Euryarchaeota'='hotpink4',
                   'Firmicutes'='brown4',
                   'Fusobacteria'='orange',
                   'Gemmatimonadetes'='darkolivegreen',
                   'Kiritimatiellaeota'='darkkhaki',
                   'Marinimicrobia_(SAR406_clade)'='darkgoldenrod4',
                   'Latescibacteria'='darkseagreen2',
                   'Lentisphaerae'='darkseagreen4',
                   'Patescibacteria'='darkturquoise',
                   'Planctomycetes'='darkslategray',
                   'Proteobacteria'='aquamarine4',
                   'Spirochaetes'='darkolivegreen3',
                   'Tenericutes'='#CA8EA7',
                   'Thaumarchaeota'='gold3',
                   'Verrucomicrobia'='darkgreen',
                   'WPS-2'='thistle2',
                   'Other'= 'black',
                   'Z-Other' = 'black')


core_physeq_genus$Phylum <- str_replace_all(core_physeq_genus$Phylum,c("Other" = "Z-Other"))
core_physeq_genus$Class <- str_replace_all(core_physeq_genus$Class,c("Other" = "Z-Other"))
save(core_physeq_genus, file = "physeq_object/core_physeq_genus.RData")

plot_phylum <- ggplot(core_physeq_genus, aes(x = diet3 , y = Abundance, fill = Phylum)) + 
  geom_bar(stat="identity", position="fill") + 
  ylab("Relative Abundance") +
  scale_y_continuous(expand = c(0,0)) + #remove the space below the 0 of the y axis in the graph
  ggtitle("Phylum Composition of Core Enteric Microbiome") +
  scale_fill_manual(values= phylum_colors) +
  theme_bw() +
  theme(axis.line = element_line(size = 0), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(), 
        axis.title = element_text(family = "serif",size = 15), 
        axis.text = element_text(size = 14), 
        axis.text.x = element_text(size = 14, family = "serif", angle = 0, vjust = 0.65), 
        axis.text.y = element_text(family = "serif", size = 13), 
        plot.title = element_text(family = "serif", size = 15), 
        legend.text = element_text(size = 10,  family = "serif"), 
        legend.title = element_text(size = 12, family = "serif"))
plot_phylum

colfunc <-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))
colors_class_rainbow <- colfunc(length(levels(as.factor(core_physeq_genus$Class))))

plot_class <- ggplot(core_physeq_genus, aes(x = diet3 , y = Abundance, fill = Class)) + 
  geom_bar(stat="identity", position="fill") + 
  ylab("Relative Abundance") +
  scale_y_continuous(expand = c(0,0)) + #remove the space below the 0 of the y axis in the graph
  ggtitle("Class Composition of Core Enteric Microbiome") +
  scale_fill_manual(values= colors_class_rainbow) +
  theme_bw() +
  theme(axis.line = element_line(size = 0), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(), 
        axis.title = element_text(family = "serif",size = 15), 
        axis.text = element_text(size = 14), 
        axis.text.x = element_text(size = 14, family = "serif", angle = 0, vjust = 0.65), 
        axis.text.y = element_text(family = "serif", size = 13), 
        plot.title = element_text(family = "serif", size = 15), 
        legend.text = element_text(size = 10,  family = "serif"), 
        legend.title = element_text(size = 12, family = "serif"))

plot_class

pdf(file = "barplots_diets.pdf", he = 7 , wi = 7)
plot_phylum
plot_class
dev.off()
