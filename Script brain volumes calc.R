## R version 4.4.2
## December 2024
## Author: Ana Cristina R. Gomes
## Modified by Maria Gustavsson, May 2025, (R version 4.4.3)

###### Brain Volumes calculation ######

#set working directory
setwd("/Users/mariagustavsson/Documents/Master thesis/data analysis")
# import data
brain_measures<-read.csv(file="brain measures_1_230_v2.csv", header=T, sep = ";")
str(brain_measures)


## compute volumes of each brain region using the elipsoid method

# Elipsoid method: Volume = ((Length × Width × Height) * (π/6)) on each brain region

# for paired regions, compute volumes for left and right and then sum volumes
# for single regions height is measured twice in each side, so an average is first computed before calculating the volume

#we take the relevant columns from the previous dataframe to be able to add the volumes in a new dataframe called brain_volumes
# dataframe to volumes
brain_volumes<-brain_measures[,colnames(brain_measures) %in% c("Dissection_ID", "Sample_ID", "Origin",
                                                               "Species", "Sex", "Weight", "SL", "Total_length")]


## TELENCEPHALON (TL) ##

# telencephalon is a paired region (i.e., it has left and right hemisphere)
# measurements of this structure were made:
#  for length and width on the top for left and right;
#  for height on either side respectively.

# first, compute volumes for the right and left TL hemisphere
#bm <- as.numeric(brain_measures) 
brain_volumes$TL_left <- (brain_measures$top_TL_left_L*brain_measures$top_TL_left_W*brain_measures$latL_TL_left_H)*(pi/6)
brain_volumes$TL_right <- (brain_measures$top_TL_right_L*brain_measures$top_TL_right_W*brain_measures$latR_TL_right_H)*(pi/6)

# check that both hemispheres are linearly related
plot(brain_volumes$TL_left, brain_volumes$TL_right)
cor(brain_volumes$TL_left, brain_volumes$TL_right, use = "complete.obs")
hist(abs(brain_volumes$TL_left - brain_volumes$TL_right))


# total volume of the telencephalon
brain_volumes$TL <- brain_volumes$TL_left+brain_volumes$TL_right

# check distribution
hist(brain_volumes$TL)




## OPTIC TECTUM (OT) ##

# optic tectum is a paired region (i.e., it has left and right hemisphere)
# measurements of this structure were made:
#  for length and width on the top for left and right;
#  for height on either side respectively.

# first, compute volumes for the right and left OT hemisphere
brain_volumes$OT_left <- (brain_measures$top_OT_left_L*brain_measures$top_OT_left_W*brain_measures$latL_OT_left_H)*(pi/6)
brain_volumes$OT_right <- (brain_measures$top_OT_right_L*brain_measures$top_OT_right_W*brain_measures$latR_OT_right_H)*(pi/6)

# check that both hemispheres are linearly related
plot(brain_volumes$OT_left, brain_volumes$OT_right)
cor(brain_volumes$OT_left, brain_volumes$OT_right, use = "complete.obs")
hist(abs(brain_volumes$OT_left - brain_volumes$OT_right))

# total volume of the optic tectum
brain_volumes$OT <- brain_volumes$OT_left+brain_volumes$OT_right

# check distribution
hist(brain_volumes$OT)




## CEREBELLUM (CB) ##

# cerebellum is a single region (i.e., it is a single structure)
# measurements of this structure were made:
#  for length and width on the top;
#  for height on either side for the same structure (i.e., need to be averaged)

# total volume of the cerebellum
# consider that you need to compute the height average of the measure took in each side
# furthermore, keep in mind that CB height needs to be /2
brain_volumes$CB <- (brain_measures$top_CB_L*brain_measures$top_CB_W*(((brain_measures$latL_CB_H.DM/2)+(brain_measures$latR_CB_H.DM/2))/2))*(pi/6)

# check distribution
hist(brain_volumes$CB)




## HYPOTHALAMUS (HP) ##

# hypothalamus is a paired region (i.e., it has left and right hemisphere)
# measurements of this structure were made:
#  for length and width on the bottom for left and right;
#  for height on either side respectively.

# first, compute volumes for the right and left HP hemisphere
brain_volumes$HP_left <- (brain_measures$bottom_HP_left_L*brain_measures$bottom_HP_left_W*brain_measures$latL_HP_left_H)*(pi/6)
brain_volumes$HP_right <- (brain_measures$bottom_HP_right_L*brain_measures$bottom_HP_right_W*brain_measures$latR_HP_right_H)*(pi/6)

# check that both hemispheres are linearly related
plot(brain_volumes$HP_left, brain_volumes$HP_right)
cor(brain_volumes$HP_left, brain_volumes$HP_right, use = "complete.obs")
hist(abs(brain_volumes$HP_left - brain_volumes$HP_right))

# total volume of the hypothalamus
brain_volumes$HP <- brain_volumes$HP_left+brain_volumes$HP_right

# check distribution
hist(brain_volumes$HP)




## OLFACTORY BULBS (OB) ##

# hypothalamus is a paired region (i.e., it has left and right hemisphere)
# measurements of this structure were made:
#  for length and width on the bottom for left and right;
#  for height on either side respectively.

# first, compute volumes for the right and left OB hemisphere
brain_volumes$OB_left <- (brain_measures$bottom_OB_left_L*brain_measures$bottom_OB_left_W*brain_measures$latL_OB_left_H)*(pi/6)
brain_volumes$OB_right <- (brain_measures$bottom_OB_right_L*brain_measures$bottom_OB_right_W*brain_measures$latR_OB_right_H)*(pi/6)

# check that both hemispheres are linearly related
plot(brain_volumes$OB_left, brain_volumes$OB_right)
cor(brain_volumes$OB_left, brain_volumes$OB_right, use = "complete.obs")
hist(abs(brain_volumes$OB_left - brain_volumes$OB_right))

# total volume of the olfactory bulbs
brain_volumes$OB <- brain_volumes$OB_left+brain_volumes$OB_right

# check distribution
hist(brain_volumes$OB)




## DORSAL MEDULA (DM) ##

# dorsal medula is a single region (i.e., it is a single structure)
# measurements of this structure were made:
#  for length and width on the bottom;
#  for height on either side for the same structure (i.e., need to be averaged)

# total volume of the dorsal medula
# consider that you need to compute the height average of the measure took in each side
brain_volumes$DM <- (brain_measures$bottom_DM_L*brain_measures$bottom_DM_W*((brain_measures$latL_DM_H+brain_measures$latR_DM_H)/2))*(pi/6)

# check distribution
hist(brain_volumes$DM)


#the results show close similarity between the regions (even if it is not 100% we assume it is correct since there could be some natural variation.)
#however, since there is no big difference we assume the data/calculations are correct so far

## TOTAL BRAIN VOLUME (brain_tot) ##

# total brain volume is the sum of all brain regions
brain_volumes$brain_tot <- brain_volumes$TL+brain_volumes$OT+brain_volumes$CB+brain_volumes$HP+brain_volumes$OB+brain_volumes$DM

# check distribution
hist(brain_volumes$brain_tot)


# subset dataset to only final volumes for analyses
brain_dataset<-brain_volumes[,colnames(brain_volumes) %in% c(
                                                             "Dissection_ID", "Sample_ID", "Origin",
                                                             "Species", "Sex", "Weight", "SL", "Total_length", "TL","OT",
                                                             "CB","HP","OB","DM","brain_tot" )]

## export dataset for analyses
write.table(brain_dataset,"brain_dataset.txt", sep="\t", row.names = FALSE)





## compute species average and sex average per species

library(dplyr)
brain_dataset_means <- brain_dataset %>% 
  group_by(Species) %>% 
  summarise(TL = mean(TL, na.rm=T),
            OT  = mean(OT, na.rm=T),
            CB  = mean(CB, na.rm=T),
            HP  = mean(HP, na.rm=T),
            OB  = mean(OB, na.rm=T),
            DM  = mean(DM, na.rm=T),
            brain_tot = mean(brain_tot, na.rm=T))
brain_dataset_means<-as.data.frame(brain_dataset_means)  
write.table(brain_dataset_means,"brain_dataset_means.txt", sep="\t", row.names = FALSE)



brain_dataset_means_sex <- brain_dataset %>% 
  group_by(Species, Sex) %>% 
  summarise(TL = mean(TL, na.rm=T),
            OT  = mean(OT, na.rm=T),
            CB  = mean(CB, na.rm=T),
            HP  = mean(HP, na.rm=T),
            OB  = mean(OB, na.rm=T),
            DM  = mean(DM, na.rm=T),
            brain_tot = mean(brain_tot, na.rm=T))
brain_dataset_means_sex<-as.data.frame(brain_dataset_means_sex)  
write.table(brain_dataset_means_sex,"brain_dataset_means_sex.txt", sep="\t", row.names = FALSE)



