#Bornean Frog Microbiome Data Analysis
#SourceTracker
#Sarah McGrath-Blaser
#January 2021


#####################Vertical Transmission###############################################
setwd("/Users/Sarah/Documents/Thesis/Sequencing_and_Analysis/SourceTracker/sourcetracker-1.0.1/")

metadata2 <- read.table('Frog_Micro_Metadata_2020.2 - Metadata_SourceTracker.txt',sep='\t',h=T,row.names=1,check=T,comment='')

vert_trans <- subset(metadata2, Env %in% c("female","male","tadpole-nest","water-terrarium",
                                           "nest-in","nest-out","leaf"))
vert_trans

otus2 <- read.table('otu_table.txt',sep='\t', header=T,row.names=1,check=T,skip=1,comment='')
otus2 <- t(as.matrix(otus2))

common.sample.ids.vert <- intersect(rownames(vert_trans), rownames(otus2))
otus_vert <- otus2[common.sample.ids.vert,]
vert_trans <- vert_trans[common.sample.ids.vert,]
vert_trans

if(length(common.sample.ids.vert) <= 1) {
  message <- paste(sprintf('Error: there are %d sample ids in common '),
                   'between the metadata file and data table')
  stop(message)
}

train.ix.vert <- which(vert_trans$SourceSink=='source')
test.ix.vert <- which(vert_trans$SourceSink=='sink')
envs_vert <- vert_trans$Env
if(is.element('Description',colnames(vert_trans))) desc <- vert_trans$Description

# load SourceTracker package
source('src/SourceTracker.r')

tune.results.vert <- tune.st(otus_vert[train.ix.vert,], envs_vert[train.ix.vert])

alpha1_vert <- tune.results.vert$best.alpha1

alpha2_vet <- tune.results.vert$best.alpha2

st_vert <- sourcetracker(otus_vert[train.ix.vert,], envs_vert[train.ix.vert])

results_vert <- predict(st_vert,otus_vert[test.ix.vert,], alpha1 = alpha1_vert,
                        alpha2 = alpha2_vet)

labels_vert <- sprintf('%s %s', envs_vert, desc)

plot(results_vert,labels_vert[test.ix.vert], type='pie',include.legend=T)

#From https://github.com/CamEJ/Plotting-Mothur-data-in-R-and-elsewhere/blob/master/SourceTracker%20plotting%20output%20-%20ggplot

PieNos = results_vert$proportions
PieNos = as.data.frame(PieNos)

#Write as a csv so you can go in and manually manipulate datasheet.
write.csv(PieNos, file = "vert_trans_sourcetracker_proportions.csv")
#Changed to long format with the headers "Tadpoles","Source_type","Proportions"

#Save as csv and read file back in to use in ggplot
PieNew <- read.csv(file = "vert_trans_sourcetracker_proportions_new.csv", header = T)
str(PieNew)

mycols <- c("darkorange2","darkgreen","gold2","dodgerblue4","darkmagenta","darkgrey","maroon")

tiff("vert_trans_pie_charts.tiff", units = "mm", width = 200, height = 100, res = 300)
#Plot pie charts
ggplot(PieNew, aes(x = "", y = Proportion, fill = Source_type)) + 
  geom_bar(stat = "identity", size = 0.5, color = "black") +
  geom_text(aes(label = "")) +
  coord_polar(theta = "y") +
  facet_wrap(~ Tadpoles, ncol = 5) +
  scale_fill_manual(values = mycols) +
  theme_void() +
  theme(strip.text.x = element_blank()) +
  labs(fill = "Source") +
  guides(fill=guide_legend(nrow=2)) +
  theme(legend.position = "bottom") +
  theme(legend.title = element_text(size=14, 
                                   face="bold")) +
  theme(legend.text = element_text(size = 12))
dev.off()

#Calculate the average percent of Unknown OTU sources for all vert_trans tadpoles
mean(PieNew$Proportion[PieNew$Source_type == "Unknown"] * 100)

##############################Environmental Transmission################################
metadata_1 <- read.table('Frog_Micro_Metadata_2020.2 - Metadata_SourceTracker_1.txt',sep='\t',h=T,row.names=1,check=T,comment='')
metadata_1_sub <- subset(metadata_1, Env=="control")
metadata_1_sub

Env_1 <- metadata_1$Env
metadata_1 <- metadata_1[ Env_1 != 'control',]
metadata_1
env_trans_1 <- subset(metadata_1, Env %in% c("female","male","tadpole-pond","water-pond",
                                          "nest-in","nest-out","leaf","water-terrarium",
                                          "tadpole-nest"))
env_trans_1

common.sample.ids.env_1 <- intersect(rownames(env_trans_1), rownames(otus2))
otus_env_1 <- otus2[common.sample.ids.env_1,]
env_trans_1 <- env_trans_1[common.sample.ids.env_1,]

if(length(common.sample.ids.env_1) <= 1) {
  message <- paste(sprintf('Error: there are %d sample ids in common '),
                   'between the metadata file and data table')
  stop(message)
}

train.ix.env_1 <- which(env_trans_1$SourceSink=='source')
test.ix.env_1 <- which(env_trans_1$SourceSink=='sink')
envs_env_1 <- env_trans_1$Env
if(is.element('Description',colnames(env_trans_1))) desc_env_1 <- env_trans_1$Description

# load SourceTracker package
source('src/SourceTracker.r')

tune.results.env_1 <- tune.st(otus_env_1[train.ix.env_1,], envs_env_1[train.ix.env_1])

alpha1_env_1 <- tune.results.env_1$best.alpha1

alpha2_env_1 <- tune.results.env_1$best.alpha2

st_env_1 <- sourcetracker(otus_env_1[train.ix.env_1,], envs_env_1[train.ix.env_1])

results_env_1 <- predict(st_env_1,otus_env_1[test.ix.env_1,], alpha1 = alpha1_env_1,
                       alpha2 = alpha2_env_1)

labels_env_1 <- sprintf('%s %s', envs_env_1, desc_env_1)

plot(results_env_1,labels_env_1[test.ix.env_1], type='pie',include.legend=T)

PieNos_env = results_env_1$proportions
PieNos_env = as.data.frame(PieNos_env)
PieNos_env
#Write as a csv so you can go in and manually manipulate datasheet.
write.csv(PieNos_env, file = "env_trans_sourcetracker_proportions.csv")
#Changed to long format with the headers "Tadpoles","Source_type","Proportions"

#Save as csv and read file back in to use in ggplot
PieNew_env <- read.csv(file = "env_trans_sourcetracker_proportions_new.csv", header = T)
str(PieNew_env)

mycols_env <- c("darkorange2","darkgreen","gold2","dodgerblue4","darkmagenta","darkseagreen2","darkgrey","deepskyblue3","maroon")

tiff("env_trans_pie_charts.tiff", units = "mm", width = 200, height = 100, res = 300)
#Plot pie charts
ggplot(PieNew_env, aes(x = "", y = Proportion, fill = Source_type)) + 
  geom_bar(stat = "identity", size = 0.5, color = "black") +
  geom_text(aes(label = "")) +
  coord_polar(theta = "y") +
  facet_wrap(~ Tadpoles, ncol = 5) +
  scale_fill_manual(values = mycols_env) +
  theme_void() +
  theme(strip.text.x = element_blank()) +
  labs(fill = "Source") +
  guides(fill=guide_legend(nrow=2)) +
  theme(legend.position = "bottom") +
  theme(legend.title = element_text(size=14, 
                                    face="bold")) +
  theme(legend.text = element_text(size = 12))
dev.off()

#Calculate the average percent of Unknown OTU sources for all env_trans tadpoles
mean(PieNew_env$Proportion[PieNew_env$Source_type == "Unknown"] * 100)