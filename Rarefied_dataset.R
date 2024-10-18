#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(tidyverse)
library(microbiome)
library(ggplot2)
library(vegan)
library(phyloseq)
library(fantaxtic)

setwd("C:/Users/sarah/OneDrive - University of Florida/Research_projects/Caecilian_project/Caecilian_bacterial_microbiome/")

##Note: when running this code the metadata import kept cutting out the header and making the first row of data the header, so followed the advice at the end of this tutorial and added a #q2:types row
#https://forum.qiime2.org/t/headers-removed-when-importing-phyloseq-object-for-decontam/8075/8
CB <- qza_to_phyloseq(features = "Outputs/table-dada2-caecilian-paired-reads.qza",tree = "Outputs/rooted-tree-caecilian-paired.qza",taxonomy = "Outputs/taxonomy-caecilian-paired.qza",metadata = "Inputs/Caecilian_Micro_Metadata_2020.2.txt")
CB
#Check metadata
sample_data(CB)

################ Basic SUMMARY of data ################
## total number of sequences per sample, sorted in order
sort(sample_sums(CB))

## lowest sample has 11,434 sequences and highest has 38,235
38235/11434
#Approximately 3 fold difference between highest and lowest library sizes

#Subset so only looking at caecilian and soil samples
CB_rare <- subset_samples(CB, Category %in% c("caecilian","soil")) 

########### Extract the otu and taxa table from only the caecilian and soil samples and give to Maria for DESeq2 analysis ##################################
# Extract abundance matrix from the phyloseq object
OTU1 = as(otu_table(CB_rare), "matrix")
# transpose if necessary
if(taxa_are_rows(CB_rare)){OTU1 <- t(OTU1)}
# Coerce to data.frame
OTUdf = as.data.frame(OTU1)
write.csv(OTUdf, file = "Outputs/Caecilian_soil_otu_table.csv")
tax_table <- tax_table(CB_rare)
write.csv(tax_table,file = "Outputs/Caecilian_soil_tax_table.csv")

#Rarefy the data to the lowest sample read value 
rarecurve(t(otu_table(CB_rare)), step = 50, cex = 0.5)
CB_rare <- rarefy_even_depth(CB_rare,rngseed = 1, sample.size = 11434, replace = F)
CB_rare
sum(taxa_sums(CB_rare))

########################## Use rarefied dataset to look at caecilian microboime only ####################################
CB_rare_caecil <- subset_samples(CB_rare, Category == "caecilian") 
CB_rare_caecil

##Bactera, archaea and unassigned  
get_taxa_unique(CB_rare_caecil, "Kingdom")
#[1] "d__Bacteria" "Unassigned"  "d__Archaea" 

#How many archaeal phyla?
caecilian_Archaea = subset_taxa(CB_rare_caecil, Kingdom == "d__Archaea")
get_taxa_unique(caecilian_Archaea, "Phylum")
#[1] "Nanoarchaeota"            "Methanobacteriota_A_1229"
#[3] "Thermoproteota"           "Halobacteriota"          
#[5] "Thermoplasmatota" 

#How many bacterial phyla?
caecilian_Bacteria = subset_taxa(CB_rare_caecil, Kingdom == "d__Bacteria")
get_taxa_unique(caecilian_Bacteria, "Phylum")
#[1] NA                          "Proteobacteria"           
#[3] "Patescibacteria"           "Chloroflexota"            
#[5] "Deinococcota"              "Dependentiae"             
#[7] "Armatimonadota"            "Fusobacteriota"           
#[9] "Verrucomicrobiota"         "Spirochaetota"            
#[11] "Cyanobacteria"             "Bacteroidota"             
#[13] "Firmicutes_D"              "Firmicutes_C"             
#[15] "Firmicutes_A"              "Firmicutes_B_370539"      
#[17] "Actinobacteriota"          "Dormibacterota"           
#[19] "Nitrospirota_A_437815"     "Acidobacteriota"          
#[21] "Bdellovibrionota_E"        "Desulfobacterota_B"       
#[23] "Desulfobacterota_G_459544" "Desulfobacterota_I"       
#[25] "Myxococcota_A_473307"      "Desulfobacterota_G_459546"
#[27] "Bdellovibrionota_C"        "Desulfobacterota_C"       
#[29] "UBA10199"                  "Planctomycetota"          
#[31] "FCPU426"                   "Eisenbacteria"            
#[33] "Firmicutes_G"              "Gemmatimonadota"          
#[35] "Firmicutes_B_370541"       "Elusimicrobiota"          
#[37] "Sumerlaeota"               "Eremiobacterota"          
#[39] "Methylomirabilota"         "Chlamydiota"              
#[41] "Myxococcota_A_437813"      "Campylobacterota"         
#[43] "Fibrobacterota" 

library(dplyr)
get_taxa_unique(CB_rare, "Phylum")
tax_table(CB_rare)[,"Phylum"] %>% unique %>% na.exclude %>% length
#Another way to see 47 total phyla in caecilian and soil samples

#### 171,510 total sequences
sum(taxa_sums(CB_rare_caecil))

#### 2214 ASVs
CB_rare_caecil

##Top 15 phyla and genera
top15ph <- sort(tapply(taxa_sums(CB_rare_caecil), tax_table(CB_rare_caecil)[, "Phylum"],sum), TRUE)[1:15]
top15g <- sort(tapply(taxa_sums(CB_rare_caecil), tax_table(CB_rare_caecil)[, "Genus"], sum), TRUE)[1:15]
top5s <- sort(tapply(taxa_sums(CB_rare_caecil), tax_table(CB_rare_caecil)[, "Species"], sum), TRUE)[1:5]
top15ph
#Proteobacteria              Bacteroidota              Firmicutes_A 
#117075                     28152                     16912 
#Firmicutes_D          Actinobacteriota           Acidobacteriota 
#12802                     12541                      5269 
#Verrucomicrobiota           Planctomycetota Desulfobacterota_G_459546 
#1548                      1226                       975 
#Myxococcota_A_473307          Campylobacterota             Chloroflexota 
#947                       811                       683 
#Cyanobacteria        Desulfobacterota_I           Patescibacteria 
#406                       367                       354 
top15g
#Acinetobacter      Pseudomonas_E_647464        Comamonas_F_589250 
#52482                      8687                      7683 
#Empedobacter_790298          Sphingobacterium          Faecalibacterium 
#7600                      3219                      3200 
#Undibacterium_570144 Stenotrophomonas_A_615274             Streptococcus 
#3026                      2710                      2158 
#Bradyrhizobium         Cryptobacteroides                Prevotella 
#2064                      1908                      1906 
#Faecousia   Sphingomicrobium_483265                    VBCG01 
#1700                      1519                      1422 
top5s

##What percentage of total sequences are these top phyla?
##Proteobacteria
Proteo = subset_taxa(CB_rare_caecil, Phylum == "Proteobacteria")
sum(taxa_sums(Proteo))
100266/171510
##Proteobacteria makes up 58% of all sequences

##Firmicutes
Firmicutes = subset_taxa(CB_rare_caecil, Phylum %in% c("Firmicutes_A","Firmicutes_D"))
sum(taxa_sums(Firmicutes))
28559/171510
##Firmicutes makes up 17% of all sequences

##Bacteroidota
Bacteroidota = subset_taxa(CB_rare_caecil, Phylum == "Bacteroidota")
sum(taxa_sums(Bacteroidota))
24166/171510
##Bacteroidota makes up 14% of all sequences

##Actinobacteriota
Actino = subset_taxa(CB_rare_caecil, Phylum == "Actinobacteriota")
sum(taxa_sums(Actino))
9675/171510
##Actinobacteria makes up 6% of all sequences

##Acidobacteriota
Acid = subset_taxa(CB_rare_caecil, Phylum == "Acidobacteriota")
sum(taxa_sums(Acid))
984/171510
##Acidobacteriota makes up 0.6% of all sequences


#########################Soil samples only####################################
SB <- subset_samples(CB_rare, Category == "soil")
sample_data(SB)
## total number of sequences per sample, sorted in order
sort(sample_sums(SB))

##Bactera, archaea and unassigned  
get_taxa_unique(SB, "Kingdom")

#How many archaeal phyla?
SB_Archae = subset_taxa(SB, Kingdom == "d__Archaea")
get_taxa_unique(SB_Archae, "Phylum")
#"Nanoarchaeota"            "Methanobacteriota_A_1229" "Thermoproteota"           "Halobacteriota"           "Thermoplasmatota" 
#5

#How many bacterial phyla?
SB_Bacteria = subset_taxa(SB, Kingdom == "d__Bacteria")
get_taxa_unique(SB_Bacteria, "Phylum")
# [1] NA                          "Proteobacteria"            "Patescibacteria"           "Chloroflexota"             "Deinococcota"             
#[6] "Dependentiae"              "Armatimonadota"            "Fusobacteriota"            "Verrucomicrobiota"         "Spirochaetota"            
#[11] "Cyanobacteria"             "Bacteroidota"              "Firmicutes_D"              "Firmicutes_C"              "Firmicutes_A"             
#[16] "Firmicutes_B_370539"       "Actinobacteriota"          "Firmicutes_B_370543"       "Dormibacterota"            "Nitrospirota_A_437815"    
#[21] "Acidobacteriota"           "Bdellovibrionota_E"        "Desulfobacterota_B"        "Desulfobacterota_G_459544" "Desulfobacterota_I"       
#[26] "Myxococcota_A_473307"      "Desulfobacterota_G_459546" "Bdellovibrionota_C"        "Desulfobacterota_C"        "UBA10199"                 
#[31] "Planctomycetota"           "FCPU426"                   "Eisenbacteria"             "Firmicutes_G"              "Gemmatimonadota"          
#[36] "Firmicutes_B_370541"       "Elusimicrobiota"           "Sumerlaeota"               "Eremiobacterota"           "Methylomirabilota"        
#[41] "Chlamydiota"               "Myxococcota_A_437813"      "Campylobacterota"          "Fibrobacterota"
#44; 43 minus NA

#### 34,302 total sequences
sum(taxa_sums(SB))
#### 2214 ASVs
SB
##Top 15 phyla and genera
top15phSB <- sort(tapply(taxa_sums(SB), tax_table(SB)[, "Phylum"],sum), TRUE)[1:15]
top15gSB <- sort(tapply(taxa_sums(SB), tax_table(SB)[, "Genus"], sum), TRUE)[1:15]
top5sSB <- sort(tapply(taxa_sums(SB), tax_table(SB)[, "Species"], sum), TRUE)[1:5]
top15phSB
#Proteobacteria           Acidobacteriota              Bacteroidota          Actinobacteriota 
#16809                      4285                      3986                      2866 
#Verrucomicrobiota              Firmicutes_A           Planctomycetota Desulfobacterota_G_459546 
#1035                       971                       955                       777 
#Myxococcota_A_473307             Chloroflexota           Patescibacteria        Desulfobacterota_B 
#582                       518                       311                       195 
#Firmicutes_C              Firmicutes_D        Desulfobacterota_I 
#188                       184                       108 
top15gSB
#Bradyrhizobium                  VBCG01             Acidoferrum             Microbacter 
#1692                    1357                    1218                     892 
#Variibacter Sphingomicrobium_483265                  WHSN01     Pseudolabrys_502538 
#892                     799                     795                     645 
#Acinetobacter              Solibacter               Gp1-AA122                  WRKU01 
#614                     428                     404                     390 
#Puia                   JC017           Clostridium_T 
#389                     387                     378 
top5sSB
#VBCG01 sp005881895 Microbacter jiangxiensis Pseudolabrys sp001426945   Acinetobacter brisouii 
#1357                      706                      645                      507 
#Solibacter usitatus_A 
#428 

##What percentage of total sequences are these top phyla?
##Proteobacteria
ProteoSB = subset_taxa(SB, Phylum == "Proteobacteria")
sum(taxa_sums(ProteoSB))
16809/34302
##Proteobacteria makes up 49% of all sequences

##Acidobacteriota
AcidSB = subset_taxa(SB, Phylum == "Acidobacteriota")
sum(taxa_sums(AcidSB))
4285/34302
##Acidobacteriota makes up 12% of all sequences

##Bacteroidota
BacteroidotaSB = subset_taxa(SB, Phylum == "Bacteroidota")
sum(taxa_sums(BacteroidotaSB))
3986/34302
##Bacteroidota makes up 11% of all sequences

##Actinobacteriota
ActinoSB = subset_taxa(SB, Phylum == "Actinobacteriota")
sum(taxa_sums(ActinoSB))
2866/34302
##Actinobacteria makes up 8% of all sequences

##Firmicutes
FirmicutesSB = subset_taxa(SB, Phylum %in% c("Firmicutes_D","Firmicutes_C","Firmicutes_A"))
sum(taxa_sums(FirmicutesSB))
1343/34302
##Firmicutes makes up 4% of all sequences


###################################################################################################

#### Alpha Diversity Comparison
CB_alpha <- microbiome::alpha(CB, index = "all")
head(CB_alpha)
CB_meta <- sample_data(CB)
CB_meta$Observed <- CB_alpha$observed
CB_meta$Shannon <- CB_alpha$diversity_shannon

CB_filtered <- filter_taxa(CB, function(x) sum(x > 2) > 2, TRUE)
CB
#No change

library(phyloseq)
library(vegan)
library(RVAideMemoire)

#subset to only look at caecilian samples
CB_caecilian <- subset_samples(CB_rare, Category == "caecilian")
CB_caecilian
sort(sample_sums(CB_caecilian))

##Are there significant differences by body location sampled?
#Jaccard
Jacc_CB <- phyloseq::distance(CB_caecilian,"jaccard")
df_CB <- as(sample_data(CB_caecilian), "data.frame")
df_CB <- na.omit(df_CB)
perm <- how(nperm = 999)
adonis2(Jacc_CB ~ Body_location, data = df_CB, permutations = perm,strata = df_CB$Caecilian_number)
#Not significantly different
#Bray-Curtis
bray_CB <- phyloseq::distance(CB_caecilian,"bray")
adonis2(bray_CB ~ Body_location, data = df_CB, permutations = perm,strata = df_CB$Caecilian_number)
#Not for bray-curtis either
library(RVAideMemoire)
packageVersion("RVAideMemoire")

##What about differences in the cutaneous microbiome between specimens?
#Jaccard
Jacc_CB <- phyloseq::distance(CB_caecilian,"jaccard")
adonis2(Jacc_CB ~ Caecilian_number, data = df_CB, permutations = perm)
#Significantly different
pairwise.perm.manova(Jacc_CB, df_CB$Caecilian_number, nperm = 999)
#Bray-Curtis
bray_CB <- phyloseq::distance(CB_caecilian,"bray")
adonis2(bray_CB ~ Caecilian_number, data = df_CB, permutations = perm)
#Significantly different
pairwise.perm.manova(bray_CB, df_CB$Caecilian_number, nperm = 999)

##Are caecilian samples significantly different from soil samples?
all_B <- subset_samples(CB_rare, Category %in% c("caecilian","soil"))
df_all_B <- as(sample_data(CB_rare), "data.frame")
df_all_B <- na.omit(df_all_B)
perm <- how(nperm = 999)
#Jaccard
Jacc_all_B <- phyloseq::distance(all_B,"jaccard")
adonis2(Jacc_all_B ~ Category, data = df_all_B, permutations = perm)
#Significantly different
pairwise.perm.manova(Jacc_CB, df_CB$Caecilian_number, nperm = 999)
#Bray-Curtis
bray_CB <- phyloseq::distance(CB_caecilian,"bray")
adonis2(bray_CB ~ Caecilian_number, data = df_CB, permutations = perm)
#Significantly different
pairwise.perm.manova(bray_CB, df_CB$Caecilian_number, nperm = 999)





bray_CB <- phyloseq::distance(CB_caecilian,"bray")
bray_ord <- ordinate(CB_caecilian,method = "PCoA","bray")
p <- plot_ordination(CB_caecilian, bray_ord, color = "Caecilian_number",shape = "Body_location") +
  geom_point(size = 5) +
  scale_shape_manual(values=c(21,22,23,24,25), name="Body location") +
  scale_color_manual(values = c("darkorchid4","goldenrod","darkgrey"), name = "Specimen number") +
  theme_q2r() +
  theme(text = element_text(size=20))

#Fixed scale_shape_manual issue using this tutorial https://github.com/joey711/phyloseq/issues/250
p
p$layers
p$layers <- p$layers[-1]
p
library(Cairo)
ggsave(filename = "/Users/sarah/OneDrive - University of Florida/Research_projects/Caecilian_project/Caecilian_bacterial_microbiome/Figures/Bray_curtis_caecilians_no_env.eps",
       plot = print(p),
       device = "eps",
       dpi = 300,
       width = 15,
       height = 12,
       unit = "cm")
ggsave(filename = "/Users/sarah/OneDrive - University of Florida/Research_projects/Caecilian_project/Caecilian_bacterial_microbiome/Figures/Bray_curtis_caecilians_no_env.tiff",
       plot = print(p),
       device = "tiff",
       dpi = 300,
       width = 15,
       height = 12,
       unit = "cm")

bray_CB <- phyloseq::distance(all_B,"bray")
bray_ord <- ordinate(all_B,method = "PCoA","bray")
p <- plot_ordination(all_B, bray_ord, color = "Caecilian_number",shape = "Body_location") +
  geom_point(size = 5) +
  scale_shape_manual(values=c(21,22,23,24,25), name="Body location") +
  scale_color_manual(values = c("darkorchid4","goldenrod","darkgrey"), name = "Specimen number") +
  theme_q2r() +
  theme(text = element_text(size=20))

#Fixed scale_shape_manual issue using this tutorial https://github.com/joey711/phyloseq/issues/250
p
p$layers
p$layers <- p$layers[-1]
p
library(Cairo)
ggsave(filename = "/Users/sarah/OneDrive - University of Florida/Research_projects/Caecilian_project/Caecilian_bacterial_microbiome/Figures/Bray_curtis_caecilians_no_env.eps",
       plot = print(p),
       device = "eps",
       dpi = 300,
       width = 15,
       height = 12,
       unit = "cm")

#####################Stacked bar plot
#install.packages("devtools")
#devtools::install_github('bbc/bbplot')
library(bbplot)

#####***Make a stacked barplot with genus and phyla data by specimen and body location

#https://grunwaldlab.github.io/analysis_of_microbiome_community_data_in_r/index.html
###Stacked bar plot by Phylum
#http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html
caecilians_phylum <- CB_caecilian %>% #Create new object
  tax_glom(taxrank = "Phylum") %>% #Condense by Phylum      
  transform_sample_counts(function(x) {x/sum(x)} ) %>% #Transform counts
  psmelt() %>% #Reformat to long form instead of wide
  filter(Abundance > 0.02) %>%  #Filter out low abundance counts
  arrange(Phylum) #Order Phyla


ggplot(caecilians_phylum, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  bbc_style() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(direction = -1) +
  scale_x_discrete(labels = NULL) +
  facet_grid(~Category,scales="free")+
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  theme(legend.position = "right", 
        legend.justification = "left",
        legend.text=element_text(size=8)) +
  guides(color = guide_legend(override.aes = list(size = 0.2)))

caecilians_genus <- CB_caecilian %>% #Create new object
  tax_glom(taxrank = "Genus") %>% #Condense by Genus      
  transform_sample_counts(function(x) {x/sum(x)} ) %>% #Transform counts
  psmelt() %>% #Reformat to long form instead of wide
  filter(Abundance > 0.02) %>%  #Filter out low abundance counts
  arrange(Genus) #Order Genera

g<- ggplot(caecilians_genus, aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  bbc_style() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(direction = -1) +
  scale_x_discrete(labels = NULL) +
  facet_grid(~Category,scales="free")+
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  theme(legend.position = "right", 
        legend.justification = "left",
        legend.text=element_text(size=8)) +
  guides(color = guide_legend(override.aes = list(size = 0.2)))
g
ggsave(filename = "/Users/sarah/OneDrive - University of Florida/Research_projects/Caecilian_project/Caecilian_bacterial_microbiome/Figures/Stacked_barplot_genera.eps",
       plot = print(g),
       device = "eps",
       dpi = 300,
       width = 15,
       height = 12,
       unit = "cm")

#Run again with x axis sample names to be sure what order the samples are in
ggplot(caecilians_genus, aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  bbc_style() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(direction = -1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(~Category,scales="free")+
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  theme(legend.position = "right", 
        legend.justification = "left",
        legend.text=element_text(size=8)) +
  guides(color = guide_legend(override.aes = list(size = 0.2)))


##Stacked barplot with soil samples included
all_B <- subset_samples(CB_rare, Category %in% c("caecilian","soil"))
sample_names(all_B)
all_genus <- all_B %>% #Create new object
  tax_glom(taxrank = "Genus") %>% #Condense by Genus      
  transform_sample_counts(function(x) {x/sum(x)} ) %>% #Transform counts
  psmelt() %>% #Reformat to long form instead of wide
  filter(Abundance > 0.02) %>%  #Filter out low abundance counts
  arrange(Genus) #Order Genera

ggplot(all_genus, aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  bbc_style() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(direction = -1) +
  scale_x_discrete(labels = NULL) +
  facet_grid(~Category,space = "free_x", scales = "free_x")+
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  theme(legend.position = "right", 
        legend.justification = "left",
        legend.text=element_text(size=8)) +
  guides(color = guide_legend(override.aes = list(size = 0.2)))

#if(!"devtools" %in% installed.packages()){
#  install.packages("devtools")
#}
#devtools::install_github("gmteunisse/fantaxtic")
#devtools::install_github("gmteunisse/ggnested")
library(fantaxtic)
library(ggnested)
top_nested <- nested_top_taxa(all_B,
                              top_tax_level = "Phylum",
                              nested_tax_level = "Genus",
                              n_top_taxa = 6, 
                              n_nested_taxa = 6)
nested_plot <- plot_nested_bar(ps_obj = top_nested$ps_obj,
                               top_level = "Phylum",
                               nested_level = "Genus")
nested_plot
ggsave(filename = "/Users/sarah/OneDrive - University of Florida/Research_projects/Caecilian_project/Caecilian_bacterial_microbiome/Figures/Stacked_barplot_nested.eps",
       plot = print(nested_plot),
       device = "eps",
       dpi = 300,
       width = 50,
       height = 25,
       unit = "cm")


jacc_CB <- phyloseq::distance(CB, "jaccard")
jacc_ord <- ordinate(CB, method = "PCoA", "jaccard")
plot_ordination(CB, jacc_ord, color = "X2") +
  geom_point(size = 5) +
  theme_bw()


jacc_CB <- phyloseq::distance(all_B, "jaccard")
jacc_ord <- ordinate(all_B, method = "PCoA", "jaccard")
plot_ordination(all_B, jacc_ord, color = "X2") +
  geom_point(size = 5) +
  theme_bw()



#Set permutations
perm <- how(nperm = 999)
#Set blocks for stratification
setBlocks(perm) <- with(df_SB, Site)
#Run permutaion test
adonis2(bray_SB ~ Day, data = df_SB, permutations = perm)
