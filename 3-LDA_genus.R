### LDA on taxa depending of the different diet
### 
### 
setwd("/Users/marie-charlottecheutin/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Carnivory")
dir.create("./LDA", recursive = T)
setwd("./LDA")

### CONDITION 1 -----
### Entre région, on teste les abondances relatives des différents GENRE bactériens entre les TROPHIC_PRED. 

### __ Europa ----
europa_core <- subset_samples(gut_core, region == "Europa")

# LEfSe preparation data 
europa.diet.genus <- tax_glom(europa_core, "Genus", NArm = TRUE) #merges species that have the same Genus + filtre Genres NA
save(europa.diet.genus, file = "europa.diet.genus.RData") 
taxo = data.frame(tax_table(europa.diet.genus))
tax_name <- paste(taxo$Domain, taxo$Phylum, taxo$Class, taxo$Order, taxo$Family, taxo$Genus,  sep="|")
tax_name2 <- c("Diet",tax_name)
save(taxo, tax_name2, file= "europa.diet.genus_taxo.RData")
europa_samp_data <- data.frame(sample_data(europa.diet.genus))
europa_diet_lefse <- data.frame(europa_samp_data[,13], stringsAsFactors = FALSE)
europa.otu = data.frame(otu_table(europa.diet.genus))

europa.lefse <- cbind(europa_diet_lefse, europa.otu)
lefse2 <- cbind(tax_name2,t(europa.lefse))
write.table(lefse2, file="europa_lefse_diet.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

## __ JDN ----
jdn_core <- subset_samples(gut_core, region == "Juan_de_nova")

# LEfSe preparation data 
jdn.diet.genus <- tax_glom(jdn_core, "Genus", NArm = TRUE) #merges species that have the same Genus + filtre Genres NA
save(jdn.diet.genus, file = "jdn.diet.genus.RData") 
taxo = data.frame(tax_table(jdn.diet.genus))
tax_name <- paste(taxo$Domain, taxo$Phylum, taxo$Class, taxo$Order, taxo$Family, taxo$Genus,  sep="|")
tax_name2 <- c("Diet",tax_name)
save(taxo, tax_name2, file= "jdn.diet.genus_taxo.RData")
jdn_samp_data <- data.frame(sample_data(jdn.diet.genus))
jdn_diet_lefse <- data.frame(jdn_samp_data[,13], stringsAsFactors = FALSE)
jdn.otu = data.frame(otu_table(jdn.diet.genus))

jdn.lefse <- cbind(jdn_diet_lefse, jdn.otu)
lefse2 <- cbind(tax_name2,t(jdn.lefse))
write.table(lefse2, file="jdn_lefse_diet.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

## __ Seychelles ----
sey_core <- subset_samples(gut_core, region == "Seychelles")

# LEfSe preparation data
sey.diet.genus <- tax_glom(sey_core, "Genus", NArm = TRUE) #merges species that have the same Genus + filtre Genres NA
save(sey.diet.genus, file = "sey.diet.genus.RData") 
taxo = data.frame(tax_table(sey.diet.genus))
tax_name <- paste(taxo$Domain, taxo$Phylum, taxo$Class, taxo$Order, taxo$Family, taxo$Genus,  sep="|")
tax_name2 <- c("Diet",tax_name)
save(taxo, tax_name2, file= "sey.diet.genus_taxo.RData")
sey_samp_data <- data.frame(sample_data(sey.diet.genus))
sey_diet_lefse <- data.frame(sey_samp_data[,13], stringsAsFactors = FALSE)
sey.otu = data.frame(otu_table(sey.diet.genus))

sey.lefse <- cbind(sey_diet_lefse, sey.otu)
lefse2 <- cbind(tax_name2,t(sey.lefse))
write.table(lefse2, file="sey_lefse_diet.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

## __ Martinique ----
yca_core <- subset_samples(gut_core, region == "Martinique")

# LEfSe preparation data
yca.diet.genus <- tax_glom(yca_core, "Genus", NArm = TRUE) #merges species that have the same Genus + filtre Genres NA
save(yca.diet.genus, file = "yca.diet.genus.RData") 
taxo = data.frame(tax_table(yca.diet.genus))
tax_name <- paste(taxo$Domain, taxo$Phylum, taxo$Class, taxo$Order, taxo$Family, taxo$Genus,  sep="|")
tax_name2 <- c("Diet",tax_name)
save(taxo, tax_name2, file= "yca.diet.genus_taxo.RData")
yca_samp_data <- data.frame(sample_data(yca.diet.genus))
yca_diet_lefse <- data.frame(yca_samp_data[,13], stringsAsFactors = FALSE)
yca.otu = data.frame(otu_table(yca.diet.genus))

yca.lefse <- cbind(yca_diet_lefse, yca.otu)
lefse2 <- cbind(tax_name2,t(yca.lefse))
write.table(lefse2, file="yca_lefse_diet.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

## __ Eparses ----
scat_core <- subset_samples(gut_core, region %in% c("Europa","Juan_de_nova"))

# LEfSe preparation data
scat.diet.genus <- tax_glom(scat_core, "Genus", NArm = TRUE) #merges species that have the same Genus + filtre Genres NA
save(scat.diet.genus, file = "scat.diet.genus.RData") 
taxo = data.frame(tax_table(scat.diet.genus))
tax_name <- paste(taxo$Domain, taxo$Phylum, taxo$Class, taxo$Order, taxo$Family, taxo$Genus,  sep="|")
tax_name2 <- c("Diet",tax_name)
save(taxo, tax_name2, file= "scat.diet.genus_taxo.RData")
scat_samp_data <- data.frame(sample_data(scat.diet.genus))
scat_diet_lefse <- data.frame(scat_samp_data[,13], stringsAsFactors = FALSE)
scat.otu = data.frame(otu_table(scat.diet.genus))

scat.lefse <- cbind(scat_diet_lefse, scat.otu)
lefse2 <- cbind(tax_name2,t(scat.lefse))
write.table(lefse2, file="scat_lefse_diet.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
save(europa_core, jdn_core, sey_core,yca_core, scat_core, file = "cores_per_regions.RData")

### CONDITION 2 -----
### On teste les abondances relatives des différentes GENRES bactériens entre les DIET3. 

### __ All -----
# LEfSe preparation data 
diet.genus <- tax_glom(gut_core, "Genus", NArm = TRUE) #merges species that have the same Genus + filtre Genres NA
save(diet.genus, file = "diet.genus.RData") 
taxo = data.frame(tax_table(diet.genus))
tax_name <- paste(taxo$Domain, taxo$Phylum, taxo$Class, taxo$Order, taxo$Family, taxo$Genus,  sep="|")
tax_name2 <- c("Diet",tax_name)
save(taxo, tax_name2, file= "diet.genus_taxo.RData")
samp_data <- data.frame(sample_data(diet.genus))
diet_lefse <- data.frame(samp_data[,12], stringsAsFactors = FALSE)
otu = data.frame(otu_table(diet.genus))

lefse <- cbind(diet_lefse, otu)
lefse2 <- cbind(tax_name2,t(lefse))
write.table(lefse2, file="lefse_diet.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

### __ Europa ----
# LEfSe preparation data 
europa_diet_lefse <- data.frame(europa_samp_data[,12], stringsAsFactors = FALSE)
europa.otu = data.frame(otu_table(europa.diet.genus))
europa.lefse <- cbind(europa_diet_lefse, europa.otu)
lefse2 <- cbind(tax_name2,t(europa.lefse))
write.table(lefse2, file="europa_lefse_diet.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

## __ JDN ----
# LEfSe preparation data 
jdn_diet_lefse <- data.frame(jdn_samp_data[,12], stringsAsFactors = FALSE)
jdn.otu = data.frame(otu_table(jdn.diet.genus))
jdn.lefse <- cbind(jdn_diet_lefse, jdn.otu)
lefse2 <- cbind(tax_name2,t(jdn.lefse))
write.table(lefse2, file="jdn_lefse_diet.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

## __ Seychelles ----
# LEfSe preparation data
sey_diet_lefse <- data.frame(sey_samp_data[,12], stringsAsFactors = FALSE)
sey.otu = data.frame(otu_table(sey.diet.genus))
sey.lefse <- cbind(sey_diet_lefse, sey.otu)
lefse2 <- cbind(tax_name2,t(sey.lefse))
write.table(lefse2, file="sey_lefse_diet.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

## __ Martinique ----
# LEfSe preparation data
yca_diet_lefse <- data.frame(yca_samp_data[,12], stringsAsFactors = FALSE)
yca.otu = data.frame(otu_table(yca.diet.genus))
yca.lefse <- cbind(yca_diet_lefse, yca.otu)
lefse2 <- cbind(tax_name2,t(yca.lefse))
write.table(lefse2, file="yca_lefse_diet.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

## __ Ocean indien ----
oi_core <- subset_samples(gut_core, region %in% c("Europa","Juan_de_nova","Seychelles"))

# LEfSe preparation data
oi.diet.genus <- tax_glom(oi_core, "Genus", NArm = TRUE) #merges species that have the same Genus + filtre Genres NA
save(oi.diet.genus, file = "oi.diet.genus.RData") 
taxo = data.frame(tax_table(oi.diet.genus))
tax_name <- paste(taxo$Domain, taxo$Phylum, taxo$Class, taxo$Order, taxo$Family, taxo$Genus,  sep="|")
tax_name2 <- c("Diet",tax_name)
save(taxo, tax_name2, file= "oi.diet.genus_taxo.RData")
oi_samp_data <- data.frame(sample_data(oi.diet.genus))
oi_diet_lefse <- data.frame(oi_samp_data[,12], stringsAsFactors = FALSE)
oi.otu = data.frame(otu_table(oi.diet.genus))

oi.lefse <- cbind(oi_diet_lefse, oi.otu)
lefse2 <- cbind(tax_name2,t(oi.lefse))
write.table(lefse2, file="oi_lefse_diet.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
save(europa_core, jdn_core, sey_core,yca_core, scat_core, oi_core, file = "cores_per_regions.RData")
