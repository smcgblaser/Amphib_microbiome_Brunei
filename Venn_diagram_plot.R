#R Code for Venn Diagram Plots
#Sarah McGrath-Blaser
#August 2020
#Run code in "Frog Microbiome Analysis.R" before running this


######################### Code for Frog Microbiome Data ###################################################
#BiocManager::install("RBGL")
#remotes::install_github("js229/Vennerable")
library(Vennerable)

##Set working directory.
setwd("/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/")

#Vertical transmission; Adult cloacas only
#I want a venn diagram of ASVs shared betweeen tadpole-nest, adults, nest-in
vert_trans <- subset_samples(fnm_cfa, Category %in% c('nest-in','adult','tadpole-nest'))

rank_names(vert_trans)
vert_trans
otu_glom_melt <- tax_glom(vert_trans, taxrank = "Species") %>% 
  psmelt() #Melt to long format filter
head(otu_glom_melt)
tail(otu_glom_melt)

#Sample data from Randall Jimenez had ASV column because he processed with DADA2 pipeline. I did Deblur and didn't have an ASV column and when processing by the most specific rank (Species) many OTUS were left out because many OTU IDs have NA as their species taxonomic description; therefore, followed the info in the tutorial below to merge taxa by OTU.
#https://joey711.github.io/phyloseq/merge.html

merged_otus <- merge_taxa(vert_trans, "OTU")
merged_otus
otu_table(merged_otus)
otu_glom_melt <- psmelt(merged_otus)

library(dplyr) #dplyr:: 

#tadpole-nest
tad_nest <- subset(otu_glom_melt, Category == "tadpole-nest") %>% 
  select(Abundance, OTU) %>% #Select only the important columns
  group_by(OTU) %>% #Combine all ASVs
  summarize(OTU_sum = sum(Abundance)) #Sum each ASV abundance count

tad_nest <- filter(tad_nest, OTU_sum > 0)#Get rid of ASVs that are 0.
tad_nest <- as.vector(tad_nest$OTU)#Create vector for ASVs IDs
tad_nest
#adult
adult_otus <- subset(otu_glom_melt, Category == "adult") %>% 
  select(Abundance, OTU) %>% #Select only the important columns
  group_by(OTU) %>% #Combine all ASVs
  summarize(OTU_sum = sum(Abundance)) #Sum each ASV abundance count

adult_otus <- filter(adult_otus, OTU_sum > 0)#Get rid of ASVs that are 0.
adult_otus <- as.vector(adult_otus$OTU)#Create vector for ASVs IDs
adult_otus
#nest-in
nest_in <- subset(otu_glom_melt, Category == "nest-in") %>% 
  select(Abundance, OTU) %>% #Select only the important columns
  group_by(OTU) %>% #Combine all ASVs
  summarize(OTU_sum = sum(Abundance)) #Sum each ASV abundance count

nest_in <- filter(nest_in, OTU_sum > 0)#Get rid of ASVs that are 0.
nest_in <- as.vector(nest_in$OTU)#Create vector for ASVs IDs
nest_in
#Create a list of the datasubsets, make sure the order in the list and names is Ok
x = list(tad_nest, adult_otus, nest_in)
names(x)[1] <- "Tadpole-nest"
names(x)[2] <- "Adults"
names(x)[3] <- "Nest-inside"

(Vstem <- Venn(x))
data <- plot(Vstem, doWeights = FALSE)

# Open a pdf file
pdf("venn_salamander.pdf", width = 8.5, height = 11) 
# 2. Create a plot
plot(Vstem, doWeights = FALSE)
# Close the pdf file
dev.off() 

################################# Vertical Transmission ####################################################
#Adults including all body locations
vert_trans_1 <- subset_samples(fnm_r_final, Category %in% c('nest-in','adult','tadpole-nest',"nest-out") & Nest_no %in% c("Nest1","Nest3","Nest4"))
sample_data(vert_trans_1)

merged_otus_1 <- merge_taxa(vert_trans_1, "OTU")
merged_otus_1
otu_table(merged_otus_1)
otu_glom_melt_1 <- psmelt(merged_otus_1)

library(dplyr) #dplyr:: 

#tadpole-nest
tad_nest <- subset(otu_glom_melt_1, Category == "tadpole-nest") %>% 
  select(Abundance, OTU) %>% #Select only the important columns
  group_by(OTU) %>% #Combine all ASVs
  summarize(OTU_sum = sum(Abundance)) #Sum each ASV abundance count

tad_nest <- filter(tad_nest, OTU_sum > 0)#Get rid of ASVs that are 0.
tad_nest <- as.vector(tad_nest$OTU)#Create vector for ASVs IDs
tad_nest
#adult
adult_otus <- subset(otu_glom_melt_1, Category == "adult") %>% 
  select(Abundance, OTU) %>% #Select only the important columns
  group_by(OTU) %>% #Combine all ASVs
  summarize(OTU_sum = sum(Abundance)) #Sum each ASV abundance count

adult_otus <- filter(adult_otus, OTU_sum > 0)#Get rid of ASVs that are 0.
adult_otus <- as.vector(adult_otus$OTU)#Create vector for ASVs IDs
adult_otus
#nest-in
nest_in <- subset(otu_glom_melt_1, Category == "nest-in") %>% 
  select(Abundance, OTU) %>% #Select only the important columns
  group_by(OTU) %>% #Combine all ASVs
  summarize(OTU_sum = sum(Abundance)) #Sum each ASV abundance count

nest_in <- filter(nest_in, OTU_sum > 0)#Get rid of ASVs that are 0.
nest_in <- as.vector(nest_in$OTU)#Create vector for ASVs IDs
nest_in
#nest-out
nest_out <- subset(otu_glom_melt_1, Category == "nest-out") %>% 
  select(Abundance, OTU) %>% #Select only the important columns
  group_by(OTU) %>% #Combine all ASVs
  summarize(OTU_sum = sum(Abundance)) #Sum each ASV abundance count

nest_out <- filter(nest_out, OTU_sum > 0)#Get rid of ASVs that are 0.
nest_out <- as.vector(nest_out$OTU)#Create vector for ASVs IDs
nest_out
#Create a list of the datasubsets, make sure the order in the list and names is Ok
x = list(tad_nest, adult_otus, nest_in, nest_out)
names(x)[1] <- "Tadpole-nest"
names(x)[2] <- "Adults"
names(x)[3] <- "Nest-inside"
names(x)[4] <- "Nest-outside"

library(Vennerable)
(Vstem <- Venn(x))
data <- plot(Vstem, doWeights = FALSE)

#################################### Environmental Transmission ############################################
env_trans <- subset_samples(fnm_r_final, Category %in% c('tadpole-pond','water-enclosure','tadpole-nest','nest-out') & Nest_no %in% c("Nest1","Nest3","Nest4"))
sample_data(env_trans)
merged_otus_envs <- merge_taxa(env_trans, "OTU")
merged_otus_envs
otu_table(merged_otus_envs)
otu_glom_melt_envs <- psmelt(merged_otus_envs)

library(dplyr) #dplyr:: 

#tadpole-nest
tad_nest <- subset(otu_glom_melt_envs, Category == "tadpole-nest") %>% 
  select(Abundance, OTU) %>% #Select only the important columns
  group_by(OTU) %>% #Combine all ASVs
  summarize(OTU_sum = sum(Abundance)) #Sum each ASV abundance count

tad_nest <- filter(tad_nest, OTU_sum > 0)#Get rid of ASVs that are 0.
tad_nest <- as.vector(tad_nest$OTU)#Create vector for ASVs IDs
tad_nest
#tadpole-pond
tad_pond <- subset(otu_glom_melt_envs, Category == "tadpole-pond") %>% 
  select(Abundance, OTU) %>% #Select only the important columns
  group_by(OTU) %>% #Combine all ASVs
  summarize(OTU_sum = sum(Abundance)) #Sum each ASV abundance count

tad_pond <- filter(tad_pond, OTU_sum > 0)#Get rid of ASVs that are 0.
tad_pond <- as.vector(tad_pond$OTU)#Create vector for ASVs IDs
tad_pond
#water_enclosure
wat_encl <- subset(otu_glom_melt_envs, Category == "water-enclosure") %>% 
  select(Abundance, OTU) %>% #Select only the important columns
  group_by(OTU) %>% #Combine all ASVs
  summarize(OTU_sum = sum(Abundance)) #Sum each ASV abundance count

wat_encl <- filter(wat_encl, OTU_sum > 0)#Get rid of ASVs that are 0.
wat_encl <- as.vector(wat_encl$OTU)#Create vector for ASVs IDs
wat_encl
#nest-out
nest_out <- subset(otu_glom_melt_envs, Category == "nest-out") %>% 
  select(Abundance, OTU) %>% #Select only the important columns
  group_by(OTU) %>% #Combine all ASVs
  summarize(OTU_sum = sum(Abundance)) #Sum each ASV abundance count

nest_out <- filter(nest_out, OTU_sum > 0)#Get rid of ASVs that are 0.
nest_out <- as.vector(nest_out$OTU)#Create vector for ASVs IDs
nest_out
#Create a list of the datasubsets, make sure the order in the list and names is Ok
x = list(tad_nest, tad_pond, wat_encl, nest_out)
names(x)[1] <- "Tadpole-nest"
names(x)[2] <- "Tadpole-pond"
names(x)[3] <- "Water"
names(x)[4] <- "Nest-outside"

library(Vennerable)
(Vstem <- Venn(x))
data <- plot(Vstem, doWeights = FALSE)



