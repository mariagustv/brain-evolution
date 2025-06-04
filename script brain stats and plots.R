## R version 4.4.2
## December 2024
## Ana Cristina R. Gomes


setwd("/Users/mariagustavsson/Documents/Master thesis/data analysis")


#### IMPORT PHYLOGENY AND INDIVIDUALS DATA ####

datanalyses_individual <- read.table("brain_dataset.txt", header = T) #import data from volumes calculation
spp_match <- read.csv(file = "abb_table.csv", header=F, sep = ";")  #import table with species names list and abbreviations

## PHYLOGENY

# import and check phylogeny
library(ape) # v5.8

#import phylogenetic tree
poecilids_tree_complete<-read.nexus("Poeciliidae.nex")
plot(poecilids_tree_complete)
poecilids_tree_complete$tip.label

#create phylogenetic tree with only the relevant species
spp_complete<-c( "Poecilia_wingei", "Poecilia_reticulata" , "Poecilia_obscura", "Micropoecilia_picta" , "Poecilia_latipinna", "Poecilia_mexicana", "Limia_perugiae", "Limia_vittata" , "Phallichthys_quadripunctatus", "Neoheterandria_elegans", "Poeciliopsis_gracilis" , "Xiphophorus_variatus", "Xiphophorus_meyeri", "Xiphophorus_nezahualcoyotl" , "Girardinus_metallicus", "Gambusia_affinis", "Gambusia_vittata" , "Xiphophorus_hellerii" , "Gambusia_sexradiata" , "Limia_nigrofasciata" , "Giardinus_metallicus" , "Limia_melanogaster" , "Poecilia_sphenops" , "Poeciliopsis_prolifica" , "Limia_versicolor" , "Heterophallus_milleri" , "Limia_dominicensis" , "Heterandria_formosa" , "Gambusia_vittata")
poecilids_tree_prunned<-drop.tip(poecilids_tree_complete, setdiff(poecilids_tree_complete$tip.label, spp_complete))
plot(poecilids_tree_prunned)
poecilids_tree_prunned

is.ultrametric(poecilids_tree_prunned)  #check phylogenetic tree if it is ultrametric, should be FALSE
# FALSE

#Add whole species name from the imported table to the new dataset
datanalyses_individual$spp_phylo<-spp_match$V2[match(datanalyses_individual$Species,spp_match$V1)]
str(datanalyses_individual)


#check data and format correctly
unique(datanalyses_individual$Species)
unique(datanalyses_individual$Sex)


#### COMPUTE AVERAGE PER SPECIES DATA ####


#### brain and morphology means
library(dplyr) 

#calculate mean for each region and create new dataset, and group by species name and sex
brain_dataset_means <- datanalyses_individual %>% 
  group_by(spp_phylo) %>% 
 # group_by(Sex) %>%   #add this code for analysis of brain diff. between sexes
    summarise(TL = mean(TL, na.rm=T),
            OT  = mean(OT, na.rm=T),
            CB  = mean(CB, na.rm=T),
            HP  = mean(HP, na.rm=T),
            OB  = mean(OB, na.rm=T),
            DM  = mean(DM, na.rm=T),
            brain_tot = mean(brain_tot, na.rm=T),
            Weight  = mean(Weight, na.rm = T))
brain_dataset_means<-as.data.frame(brain_dataset_means)  
brain_dataset_means$n_brains<-table(datanalyses_individual$spp_phylo[complete.cases(datanalyses_individual$brain_tot)])[match(brain_dataset_means$spp_phylo, names(table(datanalyses_individual$spp_phylo)))]
rownames(brain_dataset_means)<-brain_dataset_means$spp_phylo


# compute remaining variables of brain
library(phytools)
## RESIDUALS 


# to compute residuals we compute residuals from a phylogenetic regression of brain volume predicted by body mass (residuals of a regression of body mass on brain volume)
lm(log(brain_tot)~log(Weight), data=brain_dataset_means)

temp_residuals<-phyl.resid(poecilids_tree_prunned,log(setNames(brain_dataset_means$Weight,brain_dataset_means$spp_phylo)), log(setNames(brain_dataset_means$brain_tot,brain_dataset_means$spp_phylo)), method="BM")$resid
brain_dataset_means$brain_tot_resid<-temp_residuals[match(brain_dataset_means$spp_phylo, rownames(temp_residuals))]

# repeat the same for all regions
#TL
lm(log(TL)~log(Weight), data=brain_dataset_means)

temp_residuals<-phyl.resid(poecilids_tree_prunned,log(setNames(brain_dataset_means$Weight,brain_dataset_means$spp_phylo)), log(setNames(brain_dataset_means$TL,brain_dataset_means$spp_phylo)), method="BM")$resid
brain_dataset_means$brain_TL_resid<-temp_residuals[match(brain_dataset_means$spp_phylo, rownames(temp_residuals))]

#OT

lm(log(OT)~log(Weight), data=brain_dataset_means)

temp_residuals<-phyl.resid(poecilids_tree_prunned,log(setNames(brain_dataset_means$Weight,brain_dataset_means$spp_phylo)), log(setNames(brain_dataset_means$OT,brain_dataset_means$spp_phylo)), method="BM")$resid
brain_dataset_means$brain_OT_resid<-temp_residuals[match(brain_dataset_means$spp_phylo, rownames(temp_residuals))]


#CB
lm(log(CB)~log(Weight), data=brain_dataset_means)

temp_residuals<-phyl.resid(poecilids_tree_prunned,log(setNames(brain_dataset_means$Weight,brain_dataset_means$spp_phylo)), log(setNames(brain_dataset_means$CB,brain_dataset_means$spp_phylo)), method="BM")$resid
brain_dataset_means$brain_CB_resid<-temp_residuals[match(brain_dataset_means$spp_phylo, rownames(temp_residuals))]

#HP

lm(log(HP)~log(Weight), data=brain_dataset_means)

temp_residuals<-phyl.resid(poecilids_tree_prunned,log(setNames(brain_dataset_means$Weight,brain_dataset_means$spp_phylo)), log(setNames(brain_dataset_means$HP,brain_dataset_means$spp_phylo)), method="BM")$resid
brain_dataset_means$brain_HP_resid<-temp_residuals[match(brain_dataset_means$spp_phylo, rownames(temp_residuals))]


#OB

lm(log(OB)~log(Weight), data=brain_dataset_means)

temp_residuals<-phyl.resid(poecilids_tree_prunned,log(setNames(brain_dataset_means$Weight,brain_dataset_means$spp_phylo)), log(setNames(brain_dataset_means$OB,brain_dataset_means$spp_phylo)), method="BM")$resid
brain_dataset_means$brain_OB_resid<-temp_residuals[match(brain_dataset_means$spp_phylo, rownames(temp_residuals))]

#DM

lm(log(DM)~log(Weight), data=brain_dataset_means)

temp_residuals<-phyl.resid(poecilids_tree_prunned,log(setNames(brain_dataset_means$Weight,brain_dataset_means$spp_phylo)), log(setNames(brain_dataset_means$DM,brain_dataset_means$spp_phylo)), method="BM")$resid
brain_dataset_means$brain_DM_resid<-temp_residuals[match(brain_dataset_means$spp_phylo, rownames(temp_residuals))]

#### RELATION BETWEEN REGIONS AND BRAIN VOLUME  - average per species ####
brain_dataset_means #check dataset to see if the average brain volumes are calculated, and residuals are added to the dataset


## PREPARE DATASET

# log brain volume to prepare data for analysis
brain_dataset_logscale<-brain_dataset_means %>%  mutate(across(!c("spp_phylo",
                                                               "n_brains",
                                                               "brain_tot_resid","brain_TL_resid","brain_OT_resid","brain_CB_resid",
                                                               "brain_OB_resid","brain_HP_resid","brain_DM_resid"), function(x) log(x)))
brain_dataset_logscale<-brain_dataset_logscale[complete.cases(brain_dataset_logscale), ]
str(brain_dataset_logscale)

#check histogram for normal distribution
hist(brain_dataset_logscale$brain_tot)
hist(brain_dataset_logscale$TL)
hist(brain_dataset_logscale$OT)
hist(brain_dataset_logscale$CB)
hist(brain_dataset_logscale$OB)
hist(brain_dataset_logscale$HP)
hist(brain_dataset_logscale$DM)

#save dataset
write.table(brain_dataset_means,"brain_dataset_means.txt", sep="\t", row.names = FALSE)
write.table(brain_dataset_logscale, "brain_dataset_logscale.txt", sep="\t", row.names = FALSE)

#### bayesian MCMC 
library(MCMCglmm)
mcmc_priors1 <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002)))
inv.phylo<-inverseA(poecilids_tree_prunned,nodes="TIPS",scale=F)
mcmc_data <- cbind(phylo = brain_dataset_logscale$spp_phylo, brain_dataset_logscale)



# pedigree works for ultrametric trees
# for non-ultrametric trees it needs to be with the inv.phylo
# to control for phylogeny it needs to include a column called phylo with the spp ID names

# inverted matrix to account for phylogeny


# TL vs body weight

MCMCglmm_TLbmass <- MCMCglmm(TL~Weight,random=~phylo,
                                family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                                prior=mcmc_priors1,data=mcmc_data,nitt=305000,burnin=5000,
                                thin=30,verbose=FALSE)

#increased to 300 000 after complete dataset to see if it improves


summary(MCMCglmm_TLbmass)
results_pgls<-summary(MCMCglmm_TLbmass)$solutions #to save results from the model

# model assumptions check
plot(MCMCglmm_TLbmass$Sol) # fixed effects, should have no pattern (left) and widely spread (right)
effectiveSize(MCMCglmm_TLbmass$Sol)# should be between 1000 - 10.000
plot(MCMCglmm_TLbmass$VCV) # random effects variances
effectiveSize(MCMCglmm_TLbmass$VCV)



# OT vs body weight

MCMCglmm_OTbmass <- MCMCglmm(OT~Weight,random=~phylo,
                             family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                             prior=mcmc_priors1,data=mcmc_data,nitt=305000,burnin=5000,
                             thin=30,verbose=FALSE)
#increase to 300 000 after complete dataset to see if it improves


summary(MCMCglmm_OTbmass)
results_pgls<-summary(MCMCglmm_OTbmass)$solutions

# model assumptions check
plot(MCMCglmm_OTbmass$Sol) # fixed effects, should have no pattern (left) and widely spread (right)
effectiveSize(MCMCglmm_OTbmass$Sol)# should be between 1000 - 10.000
plot(MCMCglmm_OTbmass$VCV) # random effects variances
effectiveSize(MCMCglmm_OTbmass$VCV)


# CB vs body weight

MCMCglmm_CBbmass <- MCMCglmm(CB~Weight,random=~phylo,
                             family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                             prior=mcmc_priors1,data=mcmc_data,nitt=305000,burnin=5000,
                             thin=30,verbose=FALSE)
#increase to 300 000 after complete dataset to see if it improves


summary(MCMCglmm_CBbmass)
results_pgls<-summary(MCMCglmm_CBbmass)$solutions

# model assumptions check
plot(MCMCglmm_CBbmass$Sol) # fixed effects, should have no pattern (left) and widely spread (right)
effectiveSize(MCMCglmm_CBbmass$Sol)# should be between 1000 - 10.000
plot(MCMCglmm_CBbmass$VCV) # random effects variances
effectiveSize(MCMCglmm_CBbmass$VCV)


# HP vs body weight

MCMCglmm_HPbmass <- MCMCglmm(HP~Weight,random=~phylo,
                             family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                             prior=mcmc_priors1,data=mcmc_data,nitt=305000,burnin=5000,
                             thin=30,verbose=FALSE)
#increase to 300 000 after complete dataset to see if it improves


summary(MCMCglmm_HPbmass)
results_pgls<-summary(MCMCglmm_HPbmass)$solutions

# model assumptions check
plot(MCMCglmm_HPbmass$Sol) # fixed effects, should have no pattern (left) and widely spread (right)
effectiveSize(MCMCglmm_HPbmass$Sol)# should be between 1000 - 10.000
plot(MCMCglmm_HPbmass$VCV) # random effects variances
effectiveSize(MCMCglmm_HPbmass$VCV)




# OB vs body weight

MCMCglmm_OBbmass <- MCMCglmm(OB~Weight,random=~phylo,
                             family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                             prior=mcmc_priors1,data=mcmc_data,nitt=305000,burnin=5000,
                             thin=30,verbose=FALSE)
#increase to 300 000 after complete dataset to see if it improves


summary(MCMCglmm_OBbmass)
results_pgls<-summary(MCMCglmm_OBbmass)$solutions

# model assumptions check
plot(MCMCglmm_OBbmass$Sol) # fixed effects, should have no pattern (left) and widely spread (right)
effectiveSize(MCMCglmm_OBbmass$Sol)# should be between 1000 - 10.000
plot(MCMCglmm_OBbmass$VCV) # random effects variances
effectiveSize(MCMCglmm_OBbmass$VCV)


# DM vs body weight

MCMCglmm_DMbmass <- MCMCglmm(DM~Weight,random=~phylo,
                             family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                             prior=mcmc_priors1,data=mcmc_data,nitt=305000,burnin=5000,
                             thin=30,verbose=FALSE)
#increase to 300 000 after complete dataset to see if it improves


summary(MCMCglmm_DMbmass)
results_pgls<-summary(MCMCglmm_DMbmass)$solutions

# model assumptions check
plot(MCMCglmm_DMbmass$Sol) # fixed effects, should have no pattern (left) and widely spread (right)
effectiveSize(MCMCglmm_DMbmass$Sol)# should be between 1000 - 10.000
plot(MCMCglmm_DMbmass$VCV) # random effects variances
effectiveSize(MCMCglmm_DMbmass$VCV)




# tot_brain vs body weight

MCMCglmm_totb_bmass <- MCMCglmm(brain_tot~Weight,random=~phylo,
                             family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                             prior=mcmc_priors1,data=mcmc_data,nitt=305000,burnin=5000,
                             thin=30,verbose=FALSE)
#increase to 300 000 after complete dataset to see if it improves


summary(MCMCglmm_totb_bmass)
results_pgls<-summary(MCMCglmm_totb_bmass)$solutions

# model assumptions check
plot(MCMCglmm_totb_bmass$Sol) # fixed effects, should have no pattern (left) and widely spread (right)
effectiveSize(MCMCglmm_totb_bmass$Sol)# should be between 1000 - 10.000
plot(MCMCglmm_totb_bmass$VCV) # random effects variances
effectiveSize(MCMCglmm_totb_bmass$VCV)


#### non-bayesian PGLS
#test for brain volumes per species, not possible for brain volumes per sex 
#test to see if there is similar result from Bayesian (for species brain diff)

# set dataset for analyses
library(caper)
compdata_means<-comparative.data(phy=poecilids_tree_prunned,data=brain_dataset_logscale,names.col=spp_phylo,vcv=T,na.omit=T,warn.dropped=T)

PGLS_TLbmass<-pgls(TL~Weight,data=compdata_means,lambda='ML')
summary(PGLS_TLbmass)

results<-summary(PGLS_TLbmass)$coef



###### to plot #####

library(ggplot2)
library(dplyr)
library(precrec)


##TL

ggplot(brain_dataset_logscale, aes(x = Weight, y = TL, color = spp_phylo)) +
    geom_smooth(method = lm,se=T,aes(group=1),color='black')+
  #geom_smooth(method=lm, se=FALSE, fullrange=FALSE,size = 0.65, aes(group=sex, colour = sex))+
  geom_point( size = 2) +
 theme_bw() + theme(aspect.ratio=1, axis.text=element_text(size=10, colour="black"),
                     axis.title=element_text(size=15,face="bold"), axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),legend.text = element_text(size=8),legend.title = element_text(size=10), plot.title = element_text(size = 20, face = "bold"),
                     axis.ticks.length=unit(.20, "cm")) +
  guides(color = guide_legend("Species"),  shape = guide_legend("Species")) +
  labs(
    x = "Log Body mass (g)",
    y = "Log Brain volume (mm3)"
  ) + ggtitle("Telencephalon (TL)")
    


##OT

ggplot(brain_dataset_logscale, aes(x = Weight, y = OT, color = spp_phylo)) +
  geom_smooth(method = lm,se=T,aes(group=1),color='black')+
  #geom_smooth(method=lm, se=FALSE, fullrange=FALSE,size = 0.65, aes(group=sex, colour = sex))+
  geom_point( size = 2) +
  theme_bw() + theme(aspect.ratio=1, axis.text=element_text(size=10, colour="black"),
                     axis.title=element_text(size=15,face="bold"), axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),legend.text = element_text(size=8),legend.title = element_text(size=10), plot.title = element_text(size = 20, face = "bold"),
                     axis.ticks.length=unit(.20, "cm")) +
  guides(color = guide_legend("Species"),  shape = guide_legend("Species")) +
  labs(
    x = "Log Body mass (g)",
    y = "Log Brain volume (mm3)"
  )+ ggtitle("Optic Tectum (OT)")


##CB

ggplot(brain_dataset_logscale, aes(x = Weight, y = CB, color = spp_phylo)) +
  geom_smooth(method = lm,se=T,aes(group=1),color='black')+
  #geom_smooth(method=lm, se=FALSE, fullrange=FALSE,size = 0.65, aes(group=sex, colour = sex))+
  geom_point( size = 2) +
  theme_bw() + theme(aspect.ratio=1, axis.text=element_text(size=10, colour="black"),
                     axis.title=element_text(size=15,face="bold"), axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),legend.text = element_text(size=8),legend.title = element_text(size=10), plot.title = element_text(size = 20, face = "bold"),
                     axis.ticks.length=unit(.20, "cm")) +
  guides(color = guide_legend("Species"),  shape = guide_legend("Species")) +
  labs(
    x = "Log Body mass (g)",
    y = "Log Brain volume (mm3)"
  )+ ggtitle("Cerebellum (CB)")


##HP

ggplot(brain_dataset_logscale, aes(x = Weight, y = HP, color = spp_phylo)) +
  geom_smooth(method = lm,se=T,aes(group=1),color='black')+
  #geom_smooth(method=lm, se=FALSE, fullrange=FALSE,size = 0.65, aes(group=sex, colour = sex))+
  geom_point( size = 2) +
  theme_bw() + theme(aspect.ratio=1, axis.text=element_text(size=10, colour="black"),
                     axis.title=element_text(size=15,face="bold"), axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),legend.text = element_text(size=8),legend.title = element_text(size=10), plot.title = element_text(size = 20, face = "bold"),
                     axis.ticks.length=unit(.20, "cm")) +
  guides(color = guide_legend("Species"),  shape = guide_legend("Species")) +
  labs(
    x = "Log Body mass (g)",
    y = "Log Brain volume (mm3)"
  )+ ggtitle("Hypothalamus (HP)")


##OB

ggplot(brain_dataset_logscale, aes(x = Weight, y = OB, color = spp_phylo)) +
  geom_smooth(method = lm,se=T,aes(group=1),color='black')+
  #geom_smooth(method=lm, se=FALSE, fullrange=FALSE,size = 0.65, aes(group=sex, colour = sex))+
  geom_point( size = 2) +
  theme_bw() + theme(aspect.ratio=1, axis.text=element_text(size=10, colour="black"),
                     axis.title=element_text(size=15,face="bold"), axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),legend.text = element_text(size=8),legend.title = element_text(size=10), plot.title = element_text(size = 20, face = "bold"),
                     axis.ticks.length=unit(.20, "cm")) +
  guides(color = guide_legend("Species"),  shape = guide_legend("Species")) +
  labs(
    x = "Log Body mass (g)",
    y = "Log Brain volume (mm3)"
  ) + ggtitle("Olfactory bulb(OB)")


##DM

ggplot(brain_dataset_logscale, aes(x = Weight, y = DM, color = spp_phylo)) +
  geom_smooth(method = lm,se=T,aes(group=1),color='black')+
  #geom_smooth(method=lm, se=FALSE, fullrange=FALSE,size = 0.65, aes(group=sex, colour = sex))+
  geom_point( size = 2) +
  theme_bw() + theme(aspect.ratio=1, axis.text=element_text(size=10, colour="black"),
                     axis.title=element_text(size=15,face="bold"), axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),legend.text = element_text(size=8),legend.title = element_text(size=10), plot.title = element_text(size = 20, face = "bold"),
                     axis.ticks.length=unit(.20, "cm")) +
  guides(color = guide_legend("Species"),  shape = guide_legend("Species")) +
  labs(
    x = "Log Body mass (g)",
    y = "Log Brain volume (mm3)"
  ) + ggtitle("Dorsal medulla (DM)")


##total brain

ggplot(brain_dataset_logscale, aes(x = Weight, y = brain_tot, color = spp_phylo)) +
  geom_smooth(method = lm,se=T,aes(group=1),color='black')+
  #geom_smooth(method=lm, se=FALSE, fullrange=FALSE,size = 0.65, aes(group=sex, colour = sex))+
  geom_point( size = 2) +
  theme_bw() + theme(aspect.ratio=1, axis.text=element_text(size=10, colour="black"),
                     axis.title=element_text(size=15,face="bold"), axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),legend.text = element_text(size=8),legend.title = element_text(size=10), plot.title = element_text(size = 20, face = "bold"),
                     axis.ticks.length=unit(.20, "cm")) +
  guides(color = guide_legend("Species"),  shape = guide_legend("Species")) +
  labs(
    x = "Log Body mass (g)",
    y = "Log Brain volume (mm3)"
  ) + ggtitle("Total brain")




#### RELATION BETWEEN BRAIN VOLUME AND BODY MASS WITH SEX EFFECT - average per sex per species ####
#across species

# here you need to do the same process as before to compute averages, residuals, log, etc
#### brain and morphology means

library(dplyr)
brain_dataset_means_sex <- datanalyses_individual %>% 
  group_by(spp_phylo, Sex) %>% 
  summarise(TL = mean(TL, na.rm=T),
            OT  = mean(OT, na.rm=T),
            CB  = mean(CB, na.rm=T),
            HP  = mean(HP, na.rm=T),
            OB  = mean(OB, na.rm=T),
            DM  = mean(DM, na.rm=T),
            brain_tot = mean(brain_tot, na.rm=T),
            Weight  = mean(Weight, na.rm = T))
brain_dataset_means_sex<-as.data.frame(brain_dataset_means_sex)  
brain_dataset_means_sex$n_brains<-table(datanalyses_individual$spp_phylo[complete.cases(datanalyses_individual$brain_tot)])[match(brain_dataset_means_sex$spp_phylo, names(table(datanalyses_individual$spp_phylo)))]
#adding sex to the dataset

# compute remaining variables of brain
library(phytools)
## RESIDUALS 

# to compute residuals we compute residuals from a phylogenetic regression of brain volume predicted by body mass (residuals of a regression of body mass on brain volume)
lm(log(brain_tot)~log(Weight), data=brain_dataset_means_sex)

temp_residuals<-phyl.resid(poecilids_tree_prunned,log(setNames(brain_dataset_means_sex$Weight,brain_dataset_means_sex$spp_phylo)), log(setNames(brain_dataset_means_sex$brain_tot,brain_dataset_means_sex$spp_phylo)), method="BM")$resid
brain_dataset_means_sex$brain_tot_resid<-temp_residuals[match(brain_dataset_means_sex$spp_phylo, rownames(temp_residuals))]

# repeat the same for all regions
#TL
lm(log(TL)~log(Weight), data=brain_dataset_means_sex)

temp_residuals<-phyl.resid(poecilids_tree_prunned,log(setNames(brain_dataset_means_sex$Weight,brain_dataset_means_sex$spp_phylo)), log(setNames(brain_dataset_means_sex$TL,brain_dataset_means_sex$spp_phylo)), method="BM")$resid
brain_dataset_means_sex$brain_TL_resid<-temp_residuals[match(brain_dataset_means_sex$spp_phylo, rownames(temp_residuals))]

#OT

lm(log(OT)~log(Weight), data=brain_dataset_means_sex)

temp_residuals<-phyl.resid(poecilids_tree_prunned,log(setNames(brain_dataset_means_sex$Weight,brain_dataset_means_sex$spp_phylo)), log(setNames(brain_dataset_means_sex$OT,brain_dataset_means_sex$spp_phylo)), method="BM")$resid
brain_dataset_means_sex$brain_OT_resid<-temp_residuals[match(brain_dataset_means_sex$spp_phylo, rownames(temp_residuals))]


#CB
lm(log(CB)~log(Weight), data=brain_dataset_means_sex)

temp_residuals<-phyl.resid(poecilids_tree_prunned,log(setNames(brain_dataset_means_sex$Weight,brain_dataset_means_sex$spp_phylo)), log(setNames(brain_dataset_means_sex$CB,brain_dataset_means_sex$spp_phylo)), method="BM")$resid
brain_dataset_means_sex$brain_CB_resid<-temp_residuals[match(brain_dataset_means_sex$spp_phylo, rownames(temp_residuals))]

#HP

lm(log(HP)~log(Weight), data=brain_dataset_means_sex)

temp_residuals<-phyl.resid(poecilids_tree_prunned,log(setNames(brain_dataset_means_sex$Weight,brain_dataset_means_sex$spp_phylo)), log(setNames(brain_dataset_means_sex$HP,brain_dataset_means_sex$spp_phylo)), method="BM")$resid
brain_dataset_means_sex$brain_HP_resid<-temp_residuals[match(brain_dataset_means_sex$spp_phylo, rownames(temp_residuals))]


#OB

lm(log(OB)~log(Weight), data=brain_dataset_means_sex)

temp_residuals<-phyl.resid(poecilids_tree_prunned,log(setNames(brain_dataset_means_sex$Weight,brain_dataset_means_sex$spp_phylo)), log(setNames(brain_dataset_means_sex$OB,brain_dataset_means_sex$spp_phylo)), method="BM")$resid
brain_dataset_means_sex$brain_OB_resid<-temp_residuals[match(brain_dataset_means_sex$spp_phylo, rownames(temp_residuals))]

#DM

lm(log(DM)~log(Weight), data=brain_dataset_means_sex)

temp_residuals<-phyl.resid(poecilids_tree_prunned,log(setNames(brain_dataset_means_sex$Weight,brain_dataset_means_sex$spp_phylo)), log(setNames(brain_dataset_means_sex$DM,brain_dataset_means_sex$spp_phylo)), method="BM")$resid
brain_dataset_means_sex$brain_DM_resid<-temp_residuals[match(brain_dataset_means_sex$spp_phylo, rownames(temp_residuals))]

##save dataset
write.table(brain_dataset_means_sex,"brain_dataset_means_sex.txt", sep="\t", row.names = FALSE)


#############################################################################################

#### RELATION BETWEEN REGIONS AND BRAIN VOLUME  - average per species ####
brain_dataset_means_sex


## PREPARE DATASET

# log and scale regions and brain volume to analyses
brain_dataset_logscale_sex<-brain_dataset_means_sex %>%  mutate(across(!c("spp_phylo",
                                                                  "n_brains",
                                                                  "brain_tot_resid","brain_TL_resid","brain_OT_resid","brain_CB_resid",
                                                                  "brain_OB_resid","brain_HP_resid","brain_DM_resid", "Sex"), function(x) log(x)))
brain_dataset_logscale_sex<-brain_dataset_logscale_sex[complete.cases(brain_dataset_logscale_sex), ]
str(brain_dataset_logscale_sex)

hist(brain_dataset_logscale_sex$brain_tot)
hist(brain_dataset_logscale_sex$TL)
hist(brain_dataset_logscale_sex$OT)
hist(brain_dataset_logscale_sex$CB)
hist(brain_dataset_logscale_sex$OB)
hist(brain_dataset_logscale_sex$HP)
hist(brain_dataset_logscale_sex$DM)

############################################

# bayesian MCMC 
#prepare priors
mcmc_data2 <- cbind(phylo = brain_dataset_logscale_sex$spp_phylo, brain_dataset_logscale_sex)
mcmc_priors2 <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002)))

library(MCMCglmm)

# pedigree works for ultrametric trees
# for non-ultrametric trees it needs to be with the inv.phylo
# to control for phylogeny it needs to include a column called phylo with the spp ID names


# brain_tot vs body mass
# with interaction
MCMCglmm_brainmass <- MCMCglmm(brain_tot~Weight*Sex,random=~phylo+spp_phylo,
                               family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                               prior=mcmc_priors2,data=mcmc_data2,nitt=505000,burnin=5000,
                               thin=50,verbose=FALSE)
summary(MCMCglmm_brainmass)


results_brain_mass_pgls<-summary(MCMCglmm_brainmass)$solutions

# model assumptions
plot(MCMCglmm_brainmass$Sol) # fixed effects, should have no pattern (left) and widely spread (right)
effectiveSize(MCMCglmm_brainmass$Sol)# should be between 1000 - 10.000
plot(MCMCglmm_brainmass$VCV) # random effects variances
effectiveSize(MCMCglmm_brainmass$VCV)


#TL

MCMCglmm_brainmass <- MCMCglmm(TL~Weight*Sex,random=~phylo+spp_phylo,
                               family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                               prior=mcmc_priors2,data=mcmc_data2,nitt=505000,burnin=5000,
                               thin=50,verbose=FALSE)
summary(MCMCglmm_brainmass)

results_brain_mass_pgls<-summary(MCMCglmm_brainmass)$solutions

# model assumptions
plot(MCMCglmm_brainmass$Sol) # fixed effects, should have no pattern (left) and widely spread (right)
effectiveSize(MCMCglmm_brainmass$Sol)# should be between 1000 - 10.000
plot(MCMCglmm_brainmass$VCV) # random effects variances
effectiveSize(MCMCglmm_brainmass$VCV)


#OT

MCMCglmm_brainmass <- MCMCglmm(OT~Weight*Sex,random=~phylo+spp_phylo,
                               family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                               prior=mcmc_priors2,data=mcmc_data2,nitt=505000,burnin=5000,
                               thin=50,verbose=FALSE)
summary(MCMCglmm_brainmass)

results_brain_mass_pgls<-summary(MCMCglmm_brainmass)$solutions

# model assumptions
plot(MCMCglmm_brainmass$Sol) # fixed effects, should have no pattern (left) and widely spread (right)
effectiveSize(MCMCglmm_brainmass$Sol)# should be between 1000 - 10.000
plot(MCMCglmm_brainmass$VCV) # random effects variances
effectiveSize(MCMCglmm_brainmass$VCV)


#HP
MCMCglmm_brainmass <- MCMCglmm(HP~Weight*Sex,random=~phylo+spp_phylo,
                               family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                               prior=mcmc_priors2,data=mcmc_data2,nitt=505000,burnin=5000,
                               thin=50,verbose=FALSE)
summary(MCMCglmm_brainmass)

results_brain_mass_pgls<-summary(MCMCglmm_brainmass)$solutions

# model assumptions
plot(MCMCglmm_brainmass$Sol) # fixed effects, should have no pattern (left) and widely spread (right)
effectiveSize(MCMCglmm_brainmass$Sol)# should be between 1000 - 10.000
plot(MCMCglmm_brainmass$VCV) # random effects variances
effectiveSize(MCMCglmm_brainmass$VCV)


#CB
MCMCglmm_brainmass <- MCMCglmm(CB~Weight*Sex,random=~phylo+spp_phylo,
                               family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                               prior=mcmc_priors2,data=mcmc_data2,nitt=505000,burnin=5000,
                               thin=50,verbose=FALSE)
summary(MCMCglmm_brainmass)

results_brain_mass_pgls<-summary(MCMCglmm_brainmass)$solutions

# model assumptions
plot(MCMCglmm_brainmass$Sol) # fixed effects, should have no pattern (left) and widely spread (right)
effectiveSize(MCMCglmm_brainmass$Sol)# should be between 1000 - 10.000
plot(MCMCglmm_brainmass$VCV) # random effects variances
effectiveSize(MCMCglmm_brainmass$VCV)

#OB
MCMCglmm_brainmass <- MCMCglmm(OB~Weight*Sex,random=~phylo+spp_phylo,
                               family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                               prior=mcmc_priors2,data=mcmc_data2,nitt=505000,burnin=5000,
                               thin=50,verbose=FALSE)
summary(MCMCglmm_brainmass)

results_brain_mass_pgls<-summary(MCMCglmm_brainmass)$solutions

# model assumptions
plot(MCMCglmm_brainmass$Sol) # fixed effects, should have no pattern (left) and widely spread (right)
effectiveSize(MCMCglmm_brainmass$Sol)# should be between 1000 - 10.000
plot(MCMCglmm_brainmass$VCV) # random effects variances
effectiveSize(MCMCglmm_brainmass$VCV)


#DM
MCMCglmm_brainmass <- MCMCglmm(DM~Weight*Sex,random=~phylo+spp_phylo,
                               family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                               prior=mcmc_priors2,data=mcmc_data2,nitt=505000,burnin=5000,
                               thin=50,verbose=FALSE)
summary(MCMCglmm_brainmass)

results_brain_mass_pgls<-summary(MCMCglmm_brainmass)$solutions

# model assumptions
plot(MCMCglmm_brainmass$Sol) # fixed effects, should have no pattern (left) and widely spread (right)
effectiveSize(MCMCglmm_brainmass$Sol)# should be between 1000 - 10.000
plot(MCMCglmm_brainmass$VCV) # random effects variances
effectiveSize(MCMCglmm_brainmass$VCV)






###### to plot #####


library(ggplot2)
library(dplyr)


label_data <- brain_dataset_logscale_sex %>%
  group_by(Sex) %>%
  summarize(
    x = max(Weight),  # or min(Weight) for left side
    y = brain_tot_resid[which.max(Weight)]
  )

label_data$label <- ifelse(label_data$Sex == "F", "Females", "Males")

#Total brain

ggplot(brain_dataset_logscale_sex, aes(x = Weight, y = brain_tot_resid, color = spp_phylo, shape = Sex)) +
  geom_point(size = 2) +
  geom_smooth(method=lm, se=FALSE, size = 0.65, fullrange = TRUE, aes(group=Sex, linetype = Sex), show.legend = TRUE)+
  geom_text(data = label_data, aes(x = x, y = y, label = label, color = NA),  
            inherit.aes = FALSE,
            hjust = 1, vjust = -0.5, size = 3) +
  theme_bw() + theme(aspect.ratio=1, axis.text=element_text(size=10, colour="black"),
                     axis.title=element_text(size=15,face="bold"), axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),legend.text = element_text(size=8),legend.title = element_text(size=10), plot.title = element_text(size = 20, face = "bold"),
                     axis.ticks.length=unit(.20, "cm")) +

  labs(
    x = "Log Body mass (g)",
    y = "Log Brain volume (mm3)",
    color = "Species",
    shape = "Sex",
    linetype = "Sex"
  ) + ggtitle("Total brain")






# TL



ggplot(brain_dataset_logscale_sex, aes(x = Weight, y = brain_TL_resid, color = spp_phylo, shape = Sex)) +
  geom_point(size = 2) +
  geom_smooth(method=lm, se=FALSE, size = 0.65, fullrange = TRUE, aes(group=Sex, linetype = Sex), show.legend = TRUE)+
  geom_text(data = label_data, aes(x = x, y = y, label = label, color = NA),  
            inherit.aes = FALSE,
            hjust = 1, vjust = -0.5, size = 3) +
  theme_bw() + theme(aspect.ratio=1, axis.text=element_text(size=10, colour="black"),
                     axis.title=element_text(size=15,face="bold"), axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),legend.text = element_text(size=8),legend.title = element_text(size=10), plot.title = element_text(size = 20, face = "bold"),
                     axis.ticks.length=unit(.20, "cm")) +
  
  labs(
    x = "Log Body mass (g)",
    y = "Log Brain volume (mm3)",
    color = "Species",
    shape = "Sex",
    linetype = "Sex"
   ) + ggtitle("Telencephalon (TL)")




# OT



ggplot(brain_dataset_logscale_sex, aes(x = Weight, y = brain_OT_resid, color = spp_phylo, shape = Sex)) +
  geom_point(size = 2) +
  geom_smooth(method=lm, se=FALSE, size = 0.65, fullrange = TRUE, aes(group=Sex, linetype = Sex), show.legend = TRUE)+
  geom_text(data = label_data, aes(x = x, y = y, label = label, color = NA),  
            inherit.aes = FALSE,
            hjust = 1, vjust = -0.5, size = 3) +
  theme_bw() + theme(aspect.ratio=1, axis.text=element_text(size=10, colour="black"),
                     axis.title=element_text(size=15,face="bold"), axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),legend.text = element_text(size=8),legend.title = element_text(size=10), plot.title = element_text(size = 20, face = "bold"),
                     axis.ticks.length=unit(.20, "cm")) +
  
  labs(
    x = "Log Body mass (g)",
    y = "Log Brain volume (mm3)",
    color = "Species",
    shape = "Sex",
    linetype = "Sex"
  ) + ggtitle("Optic tectum(OT)")


# CB


ggplot(brain_dataset_logscale_sex, aes(x = Weight, y = brain_CB_resid, color = spp_phylo, shape = Sex)) +
  geom_point(size = 2) +
  geom_smooth(method=lm, se=FALSE, size = 0.65, fullrange = TRUE, aes(group=Sex, linetype = Sex), show.legend = TRUE)+
  geom_text(data = label_data, aes(x = x, y = y, label = label, color = NA),  
            inherit.aes = FALSE,
            hjust = 1, vjust = -0.5, size = 3) +
  theme_bw() + theme(aspect.ratio=1, axis.text=element_text(size=10, colour="black"),
                     axis.title=element_text(size=15,face="bold"), axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),legend.text = element_text(size=8),legend.title = element_text(size=10), plot.title = element_text(size = 20, face = "bold"),
                     axis.ticks.length=unit(.20, "cm")) +
  
  labs(
    x = "Log Body mass (g)",
    y = "Log Brain volume (mm3)",
    color = "Species",
    shape = "Sex",
    linetype = "Sex"
  ) + ggtitle("Cerebellum (CB)")


# HP



ggplot(brain_dataset_logscale_sex, aes(x = Weight, y = brain_HP_resid, color = spp_phylo, shape = Sex)) +
  geom_point(size = 2) +
  geom_smooth(method=lm, se=FALSE, size = 0.65, fullrange = TRUE, aes(group=Sex, linetype = Sex), show.legend = TRUE)+
  geom_text(data = label_data, aes(x = x, y = y, label = label, color = NA),  
            inherit.aes = FALSE,
            hjust = 1, vjust = -0.5, size = 3) +
  theme_bw() + theme(aspect.ratio=1, axis.text=element_text(size=10, colour="black"),
                     axis.title=element_text(size=15,face="bold"), axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),legend.text = element_text(size=8),legend.title = element_text(size=10), plot.title = element_text(size = 20, face = "bold"),
                     axis.ticks.length=unit(.20, "cm")) +
  
  labs(
    x = "Log Body mass (g)",
    y = "Log Brain volume (mm3)",
    color = "Species",
    shape = "Sex",
    linetype = "Sex"
  ) + ggtitle("Hypothalamus (HP)")

# OB



ggplot(brain_dataset_logscale_sex, aes(x = Weight, y = brain_OB_resid, color = spp_phylo, shape = Sex)) +
  geom_point(size = 2) +
  geom_smooth(method=lm, se=FALSE, size = 0.65, fullrange = TRUE, aes(group=Sex, linetype = Sex), show.legend = TRUE)+
  geom_text(data = label_data, aes(x = x, y = y, label = label, color = NA),  
            inherit.aes = FALSE,
            hjust = 1, vjust = -0.5, size = 3) +
  theme_bw() + theme(aspect.ratio=1, axis.text=element_text(size=10, colour="black"),
                     axis.title=element_text(size=15,face="bold"), axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),legend.text = element_text(size=8),legend.title = element_text(size=10), plot.title = element_text(size = 20, face = "bold"),
                     axis.ticks.length=unit(.20, "cm")) +
  
  labs(
    x = "Log Body mass (g)",
    y = "Log Brain volume (mm3)",
    color = "Species",
    shape = "Sex",
    linetype = "Sex"
  ) + ggtitle("Olfactory bulb(OB)")


# DM




ggplot(brain_dataset_logscale_sex, aes(x = Weight, y = brain_DM_resid, color = spp_phylo, shape = Sex)) +
  geom_point(size = 2) +
  geom_smooth(method=lm, se=FALSE, size = 0.65, fullrange = TRUE, aes(group=Sex, linetype = Sex), show.legend = TRUE)+
  geom_text(data = label_data, aes(x = x, y = y, label = label, color = NA),  
            inherit.aes = FALSE,
            hjust = 1, vjust = -0.5, size = 3) +
  theme_bw() + theme(aspect.ratio=1, axis.text=element_text(size=10, colour="black"),
                     axis.title=element_text(size=15,face="bold"), axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),legend.text = element_text(size=8),legend.title = element_text(size=10), plot.title = element_text(size = 20, face = "bold"),
                     axis.ticks.length=unit(.20, "cm")) +
  
  labs(
    x = "Log Body mass (g)",
    y = "Log Brain volume (mm3)",
    color = "Species",
    shape = "Sex",
    linetype = "Sex"
  ) + ggtitle("Dorsal medulla (DM)")







#Maybe not needed?
#### CORRELATION BETWEEN REGIONS ####
head(brain_dataset)

# all individuals
datanalyses_meanregion<-brain_dataset[,colnames(brain_dataset) %in% c("TL","OT","CB","HP",
                                                                      "OB","DM","spp_phylo")]
datanalyses_meanregion<-datanalyses_meanregion[complete.cases(datanalyses_meanregion), ]
str(datanalyses_meanregion)


library(phytools)
obj<-phyl.vcv(datanalyses_meanregion,vcv(poecilids_tree_prunned),1)
obj$R
## correlation between x & y
r.xy<-cov2cor(obj$R)
## t-statistic & P-value
t.xy<-r.xy*sqrt((Ntip(poecilids_tree_prunned)-2)/(1-r.xy^2))
P.xy<-2*pt(abs(t.xy),df=Ntip(poecilids_tree_prunned)-2,lower.tail=F)
P.xy


# plot the correlation results
library(corrplot) #version 0.92

corrplot(r.xy, method= "shade", type="upper",
         outline=T, pch.col="black", diag=F, tl.col = "grey40", tl.srt=60, addCoef.col = 'white')
title(main = "ALL Phylo control")


# repeat the same separating males and females
# males
datanalyses_meanregionM<-brain_dataset[brain_dataset$sex=="M", colnames(brain_dataset) %in% c("TL","OT","CB","HP",
                                                                      "OB","DM","spp_phylo")]
datanalyses_meanregionM<-datanalyses_meanregionM[complete.cases(datanalyses_meanregionM), ]
str(datanalyses_meanregionM)


objM<-phyl.vcv(datanalyses_meanregionM,vcv(poecilids_tree_prunned),1)
objM$R
## correlation between x & y
r.xy_M<-cov2cor(objM$R)
## t-statistic & P-value
t.xy_M<-r.xy_M*sqrt((Ntip(poecilids_tree_prunned)-2)/(1-r.xy_M^2))
P.xy_M<-2*pt(abs(t.xy_M),df=Ntip(poecilids_tree_prunned)-2,lower.tail=F)
P.xy_M


# plot the correlation results
corrplot(r.xy_M, method= "shade", type="upper",
         outline=T, pch.col="black", diag=F, tl.col = "grey40", tl.srt=60, addCoef.col = 'white')
title(main = "Males Phylo control")


# females
datanalyses_meanregionF<-brain_dataset[brain_dataset$sex=="F", colnames(brain_dataset) %in% c("TL","OT","CB","HP",
                                                                                              "OB","DM","spp_phylo")]
datanalyses_meanregionF<-datanalyses_meanregionF[complete.cases(datanalyses_meanregionF), ]
str(datanalyses_meanregionF)


objF<-phyl.vcv(datanalyses_meanregionF,vcv(poecilids_tree_prunned),1)
objF$R

## correlation between x & y
r.xy_F<-cov2cor(objF$R)
## t-statistic & P-value
t.xy_F<-r.xy_F*sqrt((Ntip(poecilids_tree_prunned)-2)/(1-r.xy_F^2))
P.xy_F<-2*pt(abs(t.xy_F),df=Ntip(poecilids_tree_prunned)-2,lower.tail=F)
P.xy_F


# plot the correlation results
corrplot(r.xy_F, method= "shade", type="upper",
         outline=T, pch.col="black", diag=F, tl.col = "grey40", tl.srt=60, addCoef.col = 'white')
title(main = "Females Phylo control")








