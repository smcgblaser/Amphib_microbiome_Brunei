##Frog Microbiome Analysis
##Foam nest breeding frogs (Genus:Polypedates) from Brunei, Borneo
##Analysis performed April 2020

##Intall required packages.
##Needed up update BioConductor and install phyloseq specially.https://support.bioconductor.org/p/120724/
#if(!requireNamespace("BiocManager")){
#  install.packages("BiocManager")
#}
#BiocManager::install("phyloseq")
#install.packages("ape")
#install.packages("vegan")
#install.packages("ggplot2")

##Load libraries.
library(phyloseq)
library(ape)
library(vegan)
library(ggplot2)

##Set working directory.
setwd("/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/")

########## INITIAL ANALYSES TO SHOW THAT RUN EFFECT WAS MINIMAL IN TERMS OF BETA-DIVERSITY############
### BIOM FILE ###
biom_sep_runs <- import_biom("R_files/Exported_files/exported-table-deblur-frog-forward-reads/table-with-taxonomy.biom", parseFunction = parse_taxonomy_greengenes)
sample_names(biom_sep_runs)

### MAPPING FILE ####

map_sep_runs <- import_qiime_sample_data("Metadata/Frog_Micro_Metadata_2020.2 - Sheet1 (2).tsv")
sample_names(map_sep_runs)

#####  TREE FILE  #####
tree = read.tree("R_files/Exported_files/exported-rooted-tree-frog-forward-reads/tree.nwk")

fnm <- merge_phyloseq(biom_sep_runs, tree, map_sep_runs)
fnm

### look at run effects

jaccard_run <- distance(fnm,"jaccard", binary = T)

jacc.ord_run <-ordinate(fnm, method = "PCoA", jaccard_run)

p_run <- plot_ordination(fnm, jacc.ord_run, color = "MiSeqRun")

p_run + theme_bw() + theme(text = element_text(size = 16)) + geom_point(size = 4) + 
  geom_text(aes(label = X.SampleID), size = 4 , vjust = 1.5) 

df_fnm <- as(sample_data(fnm), "data.frame")

## permanova
adonis(jaccard_run ~ MiSeqRun, data = df_fnm)

groups_run <- df_fnm[["MiSeqRun"]]
jacc_disprun <-betadisper(jaccard_run, groups_run, type=c("median"))
plot(jacc_disprun)
anova(jacc_disprun)
boxplot(jacc_disprun)

#Based on these results, then went back into QIIME2 and grouped samples by True_SampleID based on this code https://docs.qiime2.org/2020.2/plugins/available/feature-table/group/ (see also GitHub ReadMe Doc)

########################################Samples summed by True_SampleID##############################
#See QIIME2-2020.2 code for how these were grouped.
### BIOM FILE, RUNS MERGED #####
biom_grouped <- import_biom("R_files/Exported_files/exported-table-deblur-frog-forward-reads/grouped-feature-table/grouped-table-with-taxonomy.biom", parseFunction = parse_taxonomy_greengenes)
sample_names(biom_grouped)

#### MAPPING FILE ####
### We added alpha-diversity to Frog_Metadata_grouped.txt file, see at end of this script
map_grouped <- import_qiime_sample_data("Metadata/Frog_Micro_Metadata_2020.2 - Sheet3.tsv")
sample_names(map_grouped)

#####  TREE FILE  #####
tree = read.tree("R_files/Exported_files/exported-rooted-tree-frog-forward-reads/tree.nwk")

fnm_g <- merge_phyloseq(biom_grouped,tree, map_grouped)

###Need to add alpha diversity metrics to mapping file###
## extracting your otu_table table from the phyloseq object
species_site <-as(otu_table(fnm_g), "matrix")

## transpose it to get in correct format
site_species <- t(species_site)

#install.packages("picante")
library(picante)
# prune the tree to only include tips from your otu_table (might be extra)
prunedTree <- prune.sample(site_species,tree)

#Changed include.root = F to include.root = T because I kept getting a warning code
# info on what this change might mean https://www.rdocumentation.org/packages/picante/versions/1.8.1/topics/pd
#Tutorial with include.root = T http://picante.r-forge.r-project.org/picante-intro.pdf 
PD <- pd(site_species, prunedTree, include.root = T)

#need to have both alpha and df having the same column info
PD$X.SampleID <- row.names(PD)

#now merge
df_fnm_final <- merge(df_fnm, PD, by = "X.SampleID")

#Adding alpha-diversity to meta file 
# NOW, we will just read in this file from now on
write.table(df_fnm_final, "Metadata/Frog_Metadata_Final_2020.txt", sep ="\t", row.names = T)

##Re-run with new mapping file to get final phyloseq object
### BIOM FILE, RUNS MERGED #####
biom_grouped <- import_biom("R_files/Exported_files/exported-table-deblur-frog-forward-reads/grouped-feature-table/grouped-table-with-taxonomy.biom", parseFunction = parse_taxonomy_greengenes)

#### MAPPING FILE ####
### We added alpha-diversity to Frog_Metadata_grouped.txt file
map_grouped_w_PD <- import_qiime_sample_data("Metadata/Frog_Metadata_final.txt")

#####  TREE FILE  #####
tree = read.tree("R_files/Exported_files/exported-rooted-tree-frog-forward-reads/tree.nwk")

fnm_g_w_PD <- merge_phyloseq(biom_grouped,tree, map_grouped_w_PD)


############## FILTER OTUS ###################

## keep taxa that have more than two sequences in more than 1 sample 
# It has at least two sequences in at least 2 samples 

## from here: https://joey711.github.io/phyloseq/preprocess.html
## search filter_taxa
fnm_r_final = filter_taxa(fnm_g_w_PD, function(x) sum(x > 2) > 2, TRUE)
#File goes from 6.2Mb to 3.5Mb so it definitely filtered something.

## another way to do it, saying at least 2 sequences in 2% of samples (which is 2)
## This is taking the length of the otu table 
#fnm_r = filter_taxa(fnm_g, function(x) sum(x > 2) > (0.02*length(x)), TRUE)

### Reduce dataset to have only OTUs with more than 25 sequences following Bletz et al. 2017
# fnm_r = filter_taxa(fnm_g, function(x) sum(x) > 25, TRUE)

##Response from Molly Bletz about this filtering protocol:
##I chose it based on the recommendations of the deblur developer.  He said he typically removed at a 
##min of 10, but for larger datasets to bump it up to 25 or even 50 reads.  For the Mada dataset we had
##~1000 samples so I bumped it up.

################ Basic SUMMARY of data ################

## total number of sequences per sample, sorted in order
sort(sample_sums(fnm_g_w_PD))

## lowest sample has 5,244 sequences and highest has 29,267
## Approx six fold difference, which isn't too bad
29267/5244

## Based on McMurdie & Holmes 2014 we will not rarefy because of the bias in introduces and
## since our library sizes differ by six fold then according to Weiss et al. 2017 rarefying would 
## not change the false discovery rate. It only does so when library sizes are greater
## than ten fold differences

## 
get_taxa_unique(fnm_r_final, "Kingdom")

fnm_Archae = subset_taxa(fnm_r_final, Kingdom == "Archaea")

## these are your Archael phyla
get_taxa_unique(fnm_Archae, "Phylum")

library(dplyr)
## 27 bacterial and archael phyla (one NA, and 2 Archaeal phylum)
get_taxa_unique(fnm_r_final, "Phylum")
tax_table(fnm_r_final)[,"Phylum"] %>% unique %>% na.exclude %>% length

#### 1,297,428 total sequences
sum(taxa_sums(fnm_r_final))

#### 2715 sOTUs
fnm_r_final


## % Proteobacteria
Proteo = subset_taxa(fnm_r_final, Phylum == "Proteobacteria")

## total sequence count
sum(taxa_sums(fnm_r_final))

## Proteobacteria
sum(taxa_sums(Proteo))

## Proteobacteria makes up 58% of all sequences
753330/1297428

## Firmicutes
Firmicutes = subset_taxa(fnm_r_final, Phylum == "Firmicutes")
sum(taxa_sums(Firmicutes))
198561/1297428
##Firmicutes makes up 15% of all sequences

##Bacteroidetes
Bacteroidetes = subset_taxa(fnm_r_final, Phylum == "Bacteroidetes")
sum(taxa_sums(Bacteroidetes))
174121/1297428
##Bacteroidetes makes up 13% of all sequences

##% Cyanobacteria
Cyano = subset_taxa(fnm_r_final, Phylum == "Cyanobacteria")
sum(taxa_sums(Cyano))
31306/1297428
##Cyanobacteria makes up 2% of all sequences

##% Actinobacteria
Actino = subset_taxa(fnm_r_final, Phylum == "Actinobacteria")
sum(taxa_sums(Actino))
69715/1297428
##Actinobacteria makes up 5% of all sequences

##Verrucomicrobia
Verruc = subset_taxa(fnm_r_final, Phylum == "Verrucomicrobia")
sum(taxa_sums(Verruc))
26256/1297428
##Verrucomicrobia makes up 2% of all sequences


####  distribution of samples for reduced ####
#Create data frame from filtered reads for easier data manipulation.
df_fnm <- as(sample_data(fnm_r_final), "data.frame")
readsumsdf = data.frame(nreads = sort(taxa_sums(fnm_r_final), TRUE), sorted = 1:ntaxa(fnm_r_final), type = "OTUs")

readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(fnm_r_final),  TRUE), sorted = 1:nsamples(fnm_r_final), type = "Samples"))

title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")


########## look at sex effects  ###########

###### BETA-diversity  ######
adult <- subset_samples(fnm_r_final, Category == "adult")


jaccard_adult <- distance(adult,"jaccard", binary = T)
jacc.ord_adult <-ordinate(adult, method = "PCoA", jaccard_adult)

p_adult <- plot_ordination(adult, jacc.ord_adult, color = "Sex", shape = "Fert_status")
p_adult + theme_bw() + theme(text = element_text(size = 16)) + geom_point(size = 4) + 
  geom_text(aes(label = X.SampleID), size = 4 , vjust = 1.5) 

uni_adult <- distance(adult,"uunifrac")
## This is just telling us that there are more branches 5427 than the OTUs (2714)
## it's fine

uni.ord_adult <-ordinate(adult, method = "PCoA", uni_adult)

p_adult <- plot_ordination(adult, uni.ord_adult, color = "Sex", shape = "Fert_status")
p_adult + theme_bw() + theme(text = element_text(size = 16)) + geom_point(size = 4) + 
  geom_text(aes(label = X.SampleID), size = 4 , vjust = 1.5) 


bray_adult <- distance(adult,"bray")

bray.ord_adult <-ordinate(adult, method = "PCoA", bray_adult)

p_adult <- plot_ordination(adult, bray.ord_adult, color = "Sex", shape = "Fert_status")
p_adult + theme_bw() + theme(text = element_text(size = 16)) + geom_point(size = 4) + 
  geom_text(aes(label = X.SampleID), size = 4 , vjust = 1.5) 


## PERMANOVA

## make dataframe for vegan to use for analysis
df_adult <- as(sample_data(adult), "data.frame")

## we are using strata to nest body-location within adult category
## From Belden et al. 2016? Panamanian frog species host unique skin bacterial communities
## we accounted for the body location in the models using the 'strata' argument in the adonis
## function in vegan package (Oksanen et al. 2013) in R.
adonis(jaccard_adult ~ Sex, strata = df_adult$Body_location, data = df_adult)
adonis(uni_adult ~ Sex, strata = df_adult$Body_location, data = df_adult)
adonis(bray_adult ~ Sex, strata = df_adult$Body_location, data = df_adult)

## PERMDISP
groups_adult <- df_adult[["Sex"]]
jacc_dispsex <-betadisper(jaccard_adult, groups_adult, type=c("median"))
plot(jacc_dispsex)
anova(jacc_dispsex)

#### Now let's look alpha-diversity for sex
hist(df_adult$SR)
hist(df_adult$PD)

## >0.05 data is normally distributed
shapiro.test((df_adult$SR))
shapiro.test((df_adult$PD))

## > 0.05 variances are similar, meets assumption of ANOVA
# homogeniety of variance
#install.packages("car")
library(car)
leveneTest((SR) ~ Sex, data = df_adult)
leveneTest((PD) ~ Sex, data = df_adult)

aov_adult <- aov(SR ~ Sex, data = df_adult)
summary(aov_adult)

plot(SR ~ Sex, data = df_adult)

aov_adultPD <- aov(PD ~ Sex, data = df_adult)
summary(aov_adultPD)

plot(PD ~ Sex, data = df_adult)

#######################ADULTS##########################################
### Beta Diversity ###
## Looking at adults minus the unfertilized nest adults (those associated with Nest 2 removed)
adults_mun <- subset_samples(fnm_r_final, Category == "adult" & Fert_status == "fertilized")
adults_mun

jaccard_adult_mun <- distance(adults_mun,"jaccard", binary = T)
jacc.ord_adult_mun <-ordinate(adults_mun, method = "PCoA", jaccard_adult_mun)

p_adult_mun <- plot_ordination(adults_mun, jacc.ord_adult_mun, color = "Sex")
p_adult_mun + theme_bw() + theme(text = element_text(size = 16)) + geom_point(size = 4)


uni_adult_mun <- distance(adults_mun,"uunifrac")
## This is just telling us that there are more branches 5427 than the OTUs (2714)
## it's fine

uni.ord_adult_mun <-ordinate(adults_mun, method = "PCoA", uni_adult_mun)

p_adult_mun <- plot_ordination(adults_mun, uni.ord_adult_mun, color = "Sex")
p_adult_mun + theme_bw() + theme(text = element_text(size = 16)) + geom_point(size = 4)



bray_adult_mun <- distance(adults_mun,"bray")

bray.ord_adult_mun <-ordinate(adults_mun, method = "PCoA", bray_adult_mun)

p_adult_mun <- plot_ordination(adults_mun, bray.ord_adult_mun, color = "Sex")
p_adult_mun + theme_bw() + theme(text = element_text(size = 16)) + geom_point(size = 4)


## PERMANOVA

## make dataframe for vegan to use for analysis
df_adult_mun <- as(sample_data(adults_mun), "data.frame")

## we are using strata to nest body-location within adult category
## From Belden et al. 2016? Panamanian frog species host unique skin bacterial communities
## we accounted for the body location in the models using the 'strata' argument in the adonis
## function in vegan package (Oksanen et al. 2013) in R.
adonis(jaccard_adult_mun ~ Sex, strata = df_adult_mun$Body_location, data = df_adult_mun)
adonis(uni_adult_mun ~ Sex, strata = df_adult_mun$Body_location, data = df_adult_mun)
adonis(bray_adult_mun ~ Sex, strata = df_adult_mun$Body_location, data = df_adult_mun)

## PERMDISP
groups_adult_mun <- df_adult_mun[["Sex"]]
jacc_dispsex_mun <-betadisper(jaccard_adult_mun, groups_adult_mun, type=c("median"))
plot(jacc_dispsex_mun)
anova(jacc_dispsex_mun)

#### Now let's look alpha-diversity for sex
hist(df_adult_mun$SR)
hist(df_adult_mun$PD)

## >0.05 data is normally distributed
shapiro.test((df_adult_mun$SR))
shapiro.test((df_adult_mun$PD))

## > 0.05 variances are similar, meets assumption of ANOVA
# homogeniety of variance
library(car)
leveneTest((SR) ~ Sex, data = df_adult_mun)
leveneTest((PD) ~ Sex, data = df_adult_mun)

aov_adult_mun <- aov(SR ~ Sex, data = df_adult_mun)
summary(aov_adult_mun)

plot(SR ~ Sex, data = df_adult_mun)

aov_adultPD_mun <- aov(PD ~ Sex, data = df_adult_mun)
summary(aov_adultPD_mun)

plot(PD ~ Sex, data = df_adult_mun)


##################  BODY LOCATION COMPARISONS  ###########

###### BETA-diversity  ######
## This was all run above
#adult <- subset_samples(fnm_r, Category == "adult")
#jaccard_adult <- distance(adult,"jaccard", binary = T)
#jacc.ord_adult <-ordinate(adult, method = "PCoA", jaccard_adult)
#uni_adult <- distance(adult,"uunifrac")
#uni.ord_adult <-ordinate(adult, method = "PCoA", uni_adult)
#df_adult <- as(sample_data(adult), "data.frame")


p_bl <- plot_ordination(adult, jacc.ord_adult, color = "Body_location", shape = "Fert_status")
p_bl + theme_bw() + theme(text = element_text(size = 16)) + geom_point(size = 4) + 
  geom_text(aes(label = X.SampleID), size = 4 , vjust = 1.5) 

p_bl_u <- plot_ordination(adult, uni.ord_adult, color = "Body_location", shape = "Fert_status")
p_bl_u + theme_bw() + theme(text = element_text(size = 16)) + geom_point(size = 4) + 
  geom_text(aes(label = X.SampleID), size = 4 , vjust = 1.5) 

p_bl_b <- plot_ordination(adult, bray.ord_adult, color = "Body_location", shape = "Fert_status")
p_bl_b + theme_bw() + theme(text = element_text(size = 16)) + geom_point(size = 4) + 
  geom_text(aes(label = X.SampleID), size = 4 , vjust = 1.5) 

## PERMANOVA
## group medians are not different
adonis(jaccard_adult ~ Body_location, data = df_adult)
adonis(uni_adult ~ Body_location, data = df_adult)
adonis(bray_adult ~ Body_location, data = df_adult)

groups_adult <- df_adult[["Body_location"]]

## PERMDISP
jacc_dispBody_location <-betadisper(jaccard_adult, groups_adult, type=c("median"))
plot(jacc_dispBody_location)

## dispersion is not different
anova(jacc_dispBody_location)

bray_dispBody_location <-betadisper(bray_adult, groups_adult, type=c("median"))
plot(bray_dispBody_location)

#### Now let's look alpha-diversity for Body_location

## >0.05 data is normally distributed
#shapiro.test((df_adult$SR))
#shapiro.test((df_adult$PD))

## > 0.05 variances are similar, meets assumption of ANOVA
# homogeniety of variance
leveneTest((SR) ~ Body_location, data = df_adult)
leveneTest((PD) ~ Body_location, data = df_adult)

aov_bl <- aov(SR ~ Body_location, data = df_adult)
summary(aov_bl)

plot(SR ~ Body_location, data = df_adult)

aov_blPD <- aov(PD ~ Body_location, data = df_adult)
summary(aov_blPD)

plot(PD ~ Body_location, data = df_adult)

################## NEST comparison  ###################

###### BETA-diversity  ######
nest <- subset_samples(fnm_r_final, Category == "nest" & Fert_status == "fertilized")
sample_names(nest)

jaccard_nest <- distance(nest,"jaccard", binary = T)
jacc.ord_nest <-ordinate(nest, method = "PCoA", jaccard_nest)

p_nest <- plot_ordination(nest, jacc.ord_nest, color = "Nest_no", shape = "Nest_location")
p_nest + theme_bw() + theme(text = element_text(size = 16)) + geom_point(size = 4) + 
  geom_text(aes(label = X.SampleID), size = 4 , vjust = 1.5) 

p_nest + theme_bw() + theme(text = element_text(size = 16)) + geom_point(size = 4)

uni_nest <- distance(nest,"uunifrac")
## This is just telling us that there are more branches 5427 than the OTUs (2714)
## it's fine

uni.ord_nest <-ordinate(nest, method = "PCoA", uni_nest)
p_nest <- plot_ordination(nest, uni.ord_nest, color = "Fert_status", shape = "Nest_location")
p_nest + theme_bw() + theme(text = element_text(size = 16)) + geom_point(size = 4) + 
  geom_text(aes(label = X.SampleID), size = 4 , vjust = 1.5) 


bray_nest <- distance(nest,"bray")

bray.ord_nest <-ordinate(nest, method = "PCoA", bray_nest)
p_nest <- plot_ordination(nest, bray.ord_nest, color = "Fert_status", shape = "Nest_location")
p_nest + theme_bw() + theme(text = element_text(size = 16)) + geom_point(size = 4) + 
  geom_text(aes(label = X.SampleID), size = 4 , vjust = 1.5) 


## PERMANOVA

## make dataframe for vegan to use for analysis
df_nest <- as(sample_data(nest), "data.frame")


adonis(jaccard_nest ~ Nest_location,strata = df_nest$Nest_no, data = df_nest)
adonis(uni_nest ~ Nest_location, strata = df_nest$Nest_no,data = df_nest)
adonis(bray_nest ~ Nest_location, strata = df_nest$Nest_no,data = df_nest)

## PERMDISP
groups_nest <- df_nest[["Nest_location"]]
jacc_dispnest <-betadisper(jaccard_nest, groups_nest, type=c("median"))
plot(jacc_dispnest)
anova(jacc_dispnest)

### ANOSIM to look at variation within and between
nest.ano <- anosim(jaccard_nest, df_nest$Nest_location)
summary(nest.ano)
plot(nest.ano, ylab = "Dissimilarity rank")
boxplot(nest.ano$dis.rank~nest.ano$class.vec)


#### Now let's look alpha-diversity for nest_location
hist(df_nest$SR)
hist(df_nest$PD)

## >0.05 data is normally distributed
shapiro.test((df_nest$SR))

library(car)
## tried to transform in multiple ways (log, exp, asin, cube foot, squared)
## none transformed to meet normality
## will do wilcoxon signed-rank test 
## https://stats.stackexchange.com/questions/113936/what-is-the-difference-between-the-mann-whitney-and-wilcoxon-rank-sumtest
shapiro.test((df_nest$PD))

qqp((df_nest$PD), "norm")

## > 0.05 variances are similar, meets assumption of ANOVA
# homogeniety of variance
leveneTest((SR) ~ Nest_location, data = df_nest)

## We want to account for the fact that we have multiple observations
## per nest and we want to control for that between-nest variation
## https://www.r-bloggers.com/two-way-anova-with-repeated-measures/
## https://www.gribblelab.org/stats/notes/RepeatedMeasuresANOVA.pdf
aov_nest <- aov(SR ~ Nest_location + Error(Nest_no/Nest_location), data = df_nest)
summary(aov_nest)

plot(SR ~ Nest_location, data = df_nest, ylab = "OTU richness", cex.lab = 1.5)

## could not meet normality assumptions, using wilcox test
## using paired = T because of the repeated measures
wilcox.test(PD ~ Nest_location, data = df_nest, paired = T)

plot(PD ~ Nest_location, data = df_nest)


##### We now are going to move forward with cloacal microbiome  #####

## trying to filter this way, not working correctly, challenge with NAs
## going to just remove those by name

#sample_names(fnm_c)
#fnm_cf <-  subset_samples(fnm_c,Fert_status %in% c("fertilized", "NA"))

## Removing dorsal, ventral swabs from adult, unfertilized adults and nest,
## sterile-swab (no replication), randomly select one swab from inside and
## one swab from outside
#install.packages("dplyr")
library(dplyr)

'%ni%' <- Negate('%in%')
fnm_cf <- subset_samples(fnm_r_final, Body_location %ni% c("dorsal", "ventral") & 
                           X.SampleID %ni% c ('AFCN2','AFDN2', 'AFVN2',  'AMC',
                                              'AMD', 'AMV','N21', 'N22', 'N23', 'SterileSwab', 'N1in1','N1in2',
                                              'N1out2','N1out3','N3in1','N3in2','N3out2','N3out3','N4in2',
                                              'N4in3','N4out2', 'N4out3'))

sample_names(fnm_cf)

jaccard_cf <- distance(fnm_cf, "jaccard", binary = T)

jacc.ord_cf <-ordinate(fnm_cf, method = "PCoA", jaccard_cf)

p_cf <- plot_ordination(fnm_cf, jacc.ord_cf, color = "Category")
p_cf + theme_bw() + theme(text = element_text(size = 16)) + geom_point(size = 4) 
#+ geom_text(aes(label = X.SampleID), size = 4 , vjust = 1.5) 

df_cf <-as(sample_data(fnm_cf), "data.frame")

## PERMDISP
groups_cat <- df_cf[["Category"]]
jacc_dispcat <-betadisper(jaccard_cf, groups_cat, type=c("median"))
plot(jacc_dispcat)
anova(jacc_dispnest)


##########################  MAIN ANALYSES  #####################

## Now, we just want all the adult cloaca, nest, tadpoles-nest, tadpoles-pond, water-enclosure,
## leaf, and water-terrarium

## frog nest microbiome - cloaca fertilized all
fnm_cfa <- subset_samples(fnm_cf, Category %in% c("adult", "nest",  'tadpole-nest', 
                                                  'tadpole-pond', 'water-enclosure', 'leaf', 'water-terrarium'))

sample_names(fnm_cfa)
sort(sample_sums(fnm_cfa))


jaccard_cfa <- distance(fnm_cfa,"jaccard", binary = T)
jacc.ord_cfa <-ordinate(fnm_cfa, method = "PCoA", jaccard_cfa)

p_cfa <- plot_ordination(fnm_cfa, jacc.ord_cfa, color = "Category", shape = "Nest_no")
p_cfa + theme_bw() + theme(text = element_text(size = 16)) + geom_point(size = 4) + 
  geom_text(aes(label = X.SampleID), size = 4 , vjust = 1.5) 

uni_cfa <- distance(fnm_cfa,"uunifrac")
## This is just telling us that there are more branches 5427 than the OTUs (2714)
## it's fine

uni.ord_cfa <-ordinate(fnm_cfa, method = "PCoA", uni_cfa)

p_cfa <- plot_ordination(fnm_cfa, uni.ord_cfa, color = "Category", shape = "Nest_no")
p_cfa + theme_bw() + theme(text = element_text(size = 16)) + geom_point(size = 4) + 
  geom_text(aes(label = X.SampleID), size = 4 , vjust = 1.5) 


bray_cfa <- distance(fnm_cfa,"bray")

bray.ord_cfa <-ordinate(fnm_cfa, method = "PCoA", bray_cfa)

p_cfa <- plot_ordination(fnm_cfa, bray.ord_cfa, color = "Category", shape = "Nest_no")
p_cfa + theme_bw() + theme(text = element_text(size = 16)) + geom_point(size = 4) + 
  geom_text(aes(label = X.SampleID), size = 4 , vjust = 1.5) 


## PERMANOVA

## make dataframe for vegan to use for analysis
df_cfa <- as(sample_data(fnm_cfa), "data.frame")

## From Belden et al. 2016? Panamanian frog species host unique skin bacterial communities
## we accounted for the nest_no in the models using the 'strata' argument in the adonis
## function in vegan package (Oksanen et al. 2013) in R.
adonis(jaccard_cfa ~ Category, strata = df_cfa$Nest_no, data = df_cfa)
adonis(jaccard_cfa ~ Category, data = df_cfa)
#install.packages("RVAideMemoire")
library(RVAideMemoire)
pairwise.perm.manova(jaccard_cfa, df_cfa$Category, nperm = 999)

adonis(uni_cfa ~ Category, strata = df_cfa$Nest_no, data = df_cfa)

pairwise.perm.manova(uni_cfa, df_cfa$Category, nperm = 999)

adonis(bray_cfa ~ Category, strata = df_cfa$Nest_no, data = df_cfa)

pairwise.perm.manova(bray_cfa, df_cfa$Category, nperm = 999)

## PERMDISP
groups_cfa <- df_cfa[["Category"]]
jacc_disp_cfa <-betadisper(jaccard_cfa, groups_cfa, type=c("median"))
plot(jacc_disp_cfa)
anova(jacc_disp_cfa)
boxplot(jacc_disp_cfa)
TukeyHSD(jacc_disp_cfa)


#### Now let's look alpha-diversity for these 

########### summarize samples for SR and PD  #########
library(plyr)
sum_data_SR <- ddply(df_cfa, c("Category"), summarise,
                     N    = length(Category),
                     mean = mean(SR),
                     median = median(SR),
                     sd   = sd(SR),
                     se   = sd / sqrt(N))
sum_data_SR

sum_data_PD <- ddply(df_cfa, c("Category"), summarise,
                     N    = length(Category),
                     mean = mean(PD),
                     median = median(PD),
                     sd   = sd(PD),
                     se   = sd / sqrt(N))
sum_data_PD



hist(df_cfa$SR)
hist(df_cfa$PD)

## >0.05 data is normally distributed
## close enough for species richness, fits within intervals
qqp((df_cfa$SR^2), "norm")
shapiro.test((df_cfa$SR)^2)

hist((df_cfa$SR)^2)

shapiro.test((df_cfa$PD)^2)

## > 0.05 variances are similar, meets assumption of ANOVA
# homogeniety of variance
leveneTest((SR^2) ~ Category, data = df_cfa)
leveneTest((PD^2) ~ Category, data = df_cfa)

aov_cfa <- aov(SR^2 ~ Category, data = df_cfa)
summary(aov_cfa)

plot(SR^2 ~ Category, data = df_cfa)
plot(SR ~ Category, data = df_cfa)

TukeyHSD(aov_cfa)

aov_cfaPD <- aov(PD^2 ~ Category, data = df_cfa)
summary(aov_cfaPD)

plot(PD ~ Category, data = df_cfa)
plot(PD^2 ~ Category, data = df_cfa)

TukeyHSD(aov_cfaPD)

max(df_cfa$PD)
min(df_cfa$PD)
sd(df_cfa$PD)


######################### Comparison of Tadpoles before/after Environmental Interaction#################
## frog nest microbiome - tadpole comparisons (biological replicates because each tadpole swabbed is a 
## different individual)
fnm_tads <- subset_samples(fnm_cfa, Category %in% c('tadpole-nest', 'tadpole-pond'))

##only tadpoles from foam nest (FTad) and from pond (WTad)
sample_names(fnm_tads)


jaccard_tads <- distance(fnm_tads,"jaccard", binary = T)
jacc.ord_tads <-ordinate(fnm_tads, method = "PCoA", jaccard_tads)

p_tads <- plot_ordination(fnm_tads, jacc.ord_tads, color = "Nest_no", shape = "Category")
p_tads + theme_bw() + theme(text = element_text(size = 16)) + geom_point(size = 4) + 
  geom_text(aes(label = X.SampleID), size = 4 , vjust = 1.5) 

uni_tads <- distance(fnm_tads,"uunifrac")
## This is just telling us that there are more branches 5567 than the OTUs (2784)
## it's fine

uni.ord_tads <-ordinate(fnm_tads, method = "PCoA", uni_tads)

p_tads <- plot_ordination(fnm_tads, uni.ord_tads, color = "Nest_no", shape = "Category")
p_tads + theme_bw() + theme(text = element_text(size = 16)) + geom_point(size = 4) + 
  geom_text(aes(label = X.SampleID), size = 4 , vjust = 1.5) 


bray_tads <- distance(fnm_tads,"bray")

bray.ord_tads <-ordinate(fnm_tads, method = "PCoA", bray_tads)

p_tads <- plot_ordination(fnm_tads, bray.ord_tads, color = "Nest_no", shape = "Category")
p_tads + theme_bw() + theme(text = element_text(size = 16)) + geom_point(size = 4) + 
  geom_text(aes(label = X.SampleID), size = 4 , vjust = 1.5) 

## PERMANOVA

## make dataframe for vegan to use for analysis
df_tads <- as(sample_data(fnm_tads), "data.frame")


## From Belden et al. 2016? Panamanian frog species host unique skin bacterial communities
## we accounted for the nest_no in the models using the 'strata' argument in the adonis
## function in vegan package (Oksanen et al. 2013) in R.
adonis(jaccard_tads ~ Category, strata = df_tads$Nest_no, data = df_tads)
adonis(jaccard_tads ~ Category, data = df_tads)

library(RVAideMemoire)
pairwise.perm.manova(jaccard_tads, df_tads$Category, nperm = 999)

adonis(uni_tads ~ Category, strata = df_tads$Nest_no, data = df_tads)

pairwise.perm.manova(uni_tads, df_tads$Category, nperm = 999)

adonis(bray_tads ~ Category, strata = df_tads$Nest_no, data = df_tads)

pairwise.perm.manova(bray_tads, df_tads$Category, nperm = 999)

## PERMDISP
groups_tads <- df_tads[["Category"]]
jacc_disp_tads <-betadisper(jaccard_tads, groups_tads, type=c("median"))
plot(jacc_disp_tads)
anova(jacc_disp_tads)
boxplot(jacc_disp_tads)
TukeyHSD(jacc_disp_tads)

bray_disp_tads <- betadisper(bray_tads, groups_tads, type = c("median"))
plot(bray_disp_tads)

###################Creating Stacked Bar Plot Figure###############################################
library(ggplot2)

#Alpha diversity species richness by sample....
plot_richness(fnm_cfa, x="X.SampleID", measures=c("Chao1"))

#Taxanomic stacked bar plots 
#Interested in arranging via Phylum
#Look at number of features per Phylum
table(tax_table(fnm_cfa)[, "Phylum"], exclude = NULL)
#We see that there is a Phylum characterized as "NA"

#This removes Phylum NA and any Phylum characterized by ambiguous annotation.
fnm_cfa <- subset_taxa(fnm_cfa, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
table(tax_table(fnm_cfa)[, "Phylum"], exclude = NULL)
#Now we see that the NA is removed

#https://vaulot.github.io/tutorials/Phyloseq_tutorial.html
#Removed any phyla with only one feature
filterPhyla = c("SR1","TM6","GN02","OD1","Chlamydiae")
# Filter entries with unidentified Phylum.
fnm_cfa_filt_phy = subset_taxa(fnm_cfa, !Phylum %in% filterPhyla)
table(tax_table(fnm_cfa_filt_phy)[, "Phylum"], exclude = NULL)

#https://benbowlab.github.io/Hippo.html#visualizing_the_raw_sample_data 
#This works
Merged=merge_samples(fnm_cfa, "Category")
sample_data(Merged)$Cats <- factor(sample_names(Merged)) 
Merged=transform_sample_counts(Merged,function(x) 100 * x/sum(x))
plot_bar(Merged, "Cats", "Abundance", "Phylum")

############
#install.packages("devtools")
#devtools::install_github('bbc/bbplot')
library(bbplot)

#Convert phyloseq object to data frame
#https://rdrr.io/github/vmikk/metagMisc/man/phyloseq_to_df.html
#devtools::install_github("vmikk/metagMisc")
library(metagMisc)
#Create dataframe
df_1 <- phyloseq_to_df(fnm_cfa_filt_phy)
#Get rid of all columns for rank names other than Phylum
str(df_1)
write.table(df_1, "Visualizations/Filt_Phy_OTU_Table.txt", sep ="\t", row.names = T)

###Stacked bar plot by Phylum
#http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html
###Yowza###
fnm_cfa_phylum <- fnm_cfa_filt_phy %>% #Create new object
  tax_glom(taxrank = "Phylum") %>% #Condense by Phylum      
  transform_sample_counts(function(x) {x/sum(x)} ) %>% #Transform counts
  psmelt() %>% #Reformat to long form instead of wide
  filter(Abundance > 0.02) %>%  #Filter out low abundance counts        
  arrange(Phylum) #Order Phyla

colors <- 


ggplot(fnm_cfa_phylum, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "fill") +
  bbc_style() +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~Category) +
  scale_fill_viridis_d(direction = -1) +
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  theme(legend.position = "right", 
        legend.justification = "left") +
  guides(fill = guide_legend(reverse = T))


fnm_cfa_class <- fnm_cfa_filt_phy %>% #Create new object
  tax_glom(taxrank = "Class") %>% #Condense by Phylum      
  transform_sample_counts(function(x) {x/sum(x)} ) %>% #Transform counts
  psmelt() %>% #Reformat to long form instead of wide
  filter(Abundance > 0.02) %>%  #Filter out low abundance counts        
  arrange(Class) #Order Phyla

ggplot(fnm_cfa_class, aes(x = Sample, y = Abundance, fill = Class)) + 
  geom_bar(stat = "identity", position = "fill",width = 0.75) +
  bbc_style() +
  scale_fill_manual(values = phylum_colors) +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~Category) +
  scale_fill_viridis_d(direction = -1) +
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  theme(legend.position = "right", 
        legend.justification = "left") +
  guides(fill = guide_legend(reverse = T))

###########Adults - Phylum######
adult_fnm_cfa <- subset_samples(fnm_cfa_filt_phy, Category=="adult")
adult_fnm_cfa_phylum <- adult_fnm_cfa %>% #Create new object
  tax_glom(taxrank = "Phylum") %>% #Condense by Phylum      
  transform_sample_counts(function(x) {x/sum(x)} ) %>% #Transform counts
  psmelt() %>% #Reformat to long form instead of wide
  filter(Abundance > 0.02) %>%  #Filter out low abundance counts        
  arrange(Phylum) #Order Phyla

ggplot(adult_fnm_cfa_phylum, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill") +
  bbc_style() +
  scale_fill_manual(values = phylum_colors) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(direction = -1) +
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  theme(legend.position = "right", 
        legend.justification = "left") +
  guides(fill = guide_legend(reverse = T))

###Adults - Class####
adult_fnm_cfa_class <- adult_fnm_cfa %>% #Create new object
  tax_glom(taxrank = "Class") %>% #Condense by Phylum      
  transform_sample_counts(function(x) {x/sum(x)} ) %>% #Transform counts
  psmelt() %>% #Reformat to long form instead of wide
  filter(Abundance > 0.02) %>%  #Filter out low abundance counts        
  arrange(Class) #Order Phyla

ggplot(adult_fnm_cfa_class, aes(x = Sample, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity", position = "fill") +
  bbc_style() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(direction = -1) +
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  theme(legend.position = "right", 
        legend.justification = "left") +
  guides(fill = guide_legend(reverse = T)) +
  scale_x_discrete(labels=c("female","female","female","male","male","male")) +
  theme(axis.text.x=element_text(angle=45, hjust = 0.5))

###Nest by Location - Class###
nest_fnm_cfa <- subset_samples(fnm_cfa_filt_phy, Category=="nest")
#https://github.com/joey711/phyloseq/issues/293
sample_variables(nest_fnm_cfa)
variable1 = as.character(get_variable(nest_fnm_cfa, "Nest_location"))
variable2 = as.character(get_variable(nest_fnm_cfa, "Category"))
sample_data(nest_fnm_cfa)$NestxLocation <- mapply(paste0, variable1, variable2,
                                                  collapse = "_")
merge_samples(nest_fnm_cfa, "NestxLocation")


nest_fnm_cfa_class <- nest_fnm_cfa %>% #Create new object
  tax_glom(taxrank = "Class") %>% #Condense by Class      
  transform_sample_counts(function(x) {x/sum(x)} ) %>% #Transform counts
  psmelt() %>% #Reformat to long form instead of wide
  filter(Abundance > 0.02) %>%  #Filter out low abundance counts        
  arrange(Class) #Order Class

x.ordered <- factor(nest_fnm_cfa_class$NestxLocation, levels=c("insidenest", "insidenest", 
                                            "insidenest", "outsidenest",
                                            "outsidenest", "outsidenest"))

ggplot(nest_fnm_cfa_class, aes(x = Sample, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity", position = "fill") +
  bbc_style() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(direction = -1) +
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  theme(legend.position = "right", 
        legend.justification = "left") +
  guides(fill = guide_legend(reverse = T)) +
  scale_x_discrete(labels=c("inside","outside","inside","outside","inside"
                            ,"outside")) +
  theme(axis.text.x=element_text(angle=45, hjust = 0.5))





#Network plots
plot_net(fnm_cfa, "bray", color = "Nest_no",laymeth='auto')

#Other good tutorials
#https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html#supervised_learning

