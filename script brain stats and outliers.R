## R version 4.4.2
## December 2024
## Author: Ana Cristina R. Gomes
##Modified by Maria Gustavsson, May 2025, (R version 4.4.3)


#import the necessary packages
library(dplyr)
library(tidyverse)
library(ggplot2)


#import the data
brain_dataset<-read.table(file="brain_dataset.txt", sep="\t", header=T)


#### CHECK OUTLIERS FOR EACH REGION WITHIN SPECIES ####

species_list<-unique(brain_dataset$Species)
# species
# [1] "PRE" "MPI" "PSA" "XNE" "GVI" "XVA" "PQA" "HFO" "NEL" "PMX" "LVI" "POB" "XMY" "PLA" "LPE" "PWI" "PGR" "GME" "GAF"

#brain_dataset<-brain_dataset[brain_dataset$species!="HFO",]

library(tidyr)


# Reshape the dataset to have one column with all of the brain regions
reshape_data <- brain_dataset %>%
  pivot_longer(
    cols =  c(TL, OT, CB, HP, DM, OB, DM, brain_tot), 
    names_to = "brain_region",
    values_to = "brain_volume"
  )


#Plot outliers for all species by each brain region
for(r in unique(reshape_data$brain_region)) {
  print(
    ggplot(subset(reshape_data, brain_region == r), aes(x = Species, y = brain_volume, fill = Species)) +
      geom_boxplot(outlier.colour = "red", outlier.shape = 16) +
      theme_minimal() +
      labs(title = paste("Outliers by Species for", r))
  )
}

#################################################################################
### test for outliers each species separately


## PRE
PRE_brain<-brain_dataset[brain_dataset$Species=="PRE",]

boxplot(PRE_brain$TL)
boxplot.stats(PRE_brain$TL)$out

boxplot(PRE_brain$OT)
boxplot.stats(PRE_brain$OT)$out #One outlier in OT
outliers_data<-rbind(outliers_data, unlist(c(PRE_brain[which(PRE_brain$OT %in% c(boxplot.stats(PRE_brain$OT)$out)), ], "OT")))


boxplot(PRE_brain$CB)
boxplot.stats(PRE_brain$CB)$out

boxplot(PRE_brain$HP) #One outlier in HP
boxplot.stats(PRE_brain$HP)$out
outliers_data<-cbind(PRE_brain[which(PRE_brain$HP %in% c(boxplot.stats(PRE_brain$HP)$out)), ], "HP")
colnames(outliers_data)[length(outliers_data)]<-"region"

boxplot(PRE_brain$OB)
boxplot.stats(PRE_brain$OB)$out
outliers_data<-rbind(outliers_data, unlist(c(PRE_brain[which(PRE_brain$OB %in% c(boxplot.stats(PRE_brain$OB)$out)), ], "OB")))

boxplot(PRE_brain$DM)
boxplot.stats(PRE_brain$DM)$out





## MPI
MPI_brain<-brain_dataset[brain_dataset$Species=="MPI",]

boxplot(MPI_brain$TL)
boxplot.stats(MPI_brain$TL)$out
outliers_data<-rbind(outliers_data, unlist(c(MPI_brain[which(MPI_brain$TL %in% c(boxplot.stats(MPI_brain$TL)$out)), ], "TL")))

boxplot(MPI_brain$OT)
boxplot.stats(MPI_brain$OT)$out

boxplot(MPI_brain$CB)
boxplot.stats(MPI_brain$CB)$out

boxplot(MPI_brain$HP)
boxplot.stats(MPI_brain$HP)$out
outliers_data<-rbind(outliers_data, unlist(c(MPI_brain[which(MPI_brain$HP %in% c(boxplot.stats(MPI_brain$HP)$out)), ], "HP")))

boxplot(MPI_brain$OB)
boxplot.stats(MPI_brain$OB)$out

boxplot(MPI_brain$DM)
boxplot.stats(MPI_brain$DM)$out
outliers_data<-rbind(outliers_data, unlist(c(MPI_brain[which(MPI_brain$DM %in% c(boxplot.stats(MPI_brain$DM)$out)), ], "DM")))



#PSP

PSP_brain<-brain_dataset[brain_dataset$Species=="PSP",]

boxplot(PSP_brain$TL)
boxplot.stats(PSP_brain$TL)$out

boxplot(PSP_brain$OT)
boxplot.stats(PSP_brain$OT)$out

boxplot(PSP_brain$CB)
boxplot.stats(PSP_brain$CB)$out

boxplot(PSP_brain$HP)
boxplot.stats(PSP_brain$HP)$out
outliers_data<-rbind(outliers_data, unlist(c(PSP_brain[which(PSP_brain$HP %in% c(boxplot.stats(PSP_brain$HP)$out)), ], "HP")))

boxplot(PSP_brain$OB)
boxplot.stats(PSP_brain$OB)$out
outliers_data<-rbind(outliers_data, unlist(c(PSP_brain[which(PSP_brain$OB %in% c(boxplot.stats(PSP_brain$OB)$out)), ], "OB")))


boxplot(PSP_brain$DM)
boxplot.stats(PSP_brain$DM)$out


## XNE
XNE_brain<-brain_dataset[brain_dataset$Species=="XNE",]

boxplot(XNE_brain$TL)
boxplot.stats(XNE_brain$TL)$out
outliers_data<-rbind(outliers_data, unlist(c(XNE_brain[which(XNE_brain$TL %in% c(boxplot.stats(XNE_brain$TL)$out)), ], "TL")))

boxplot(XNE_brain$OT)
boxplot.stats(XNE_brain$OT)$out
outliers_data<-rbind(outliers_data, unlist(c(XNE_brain[which(XNE_brain$OT %in% c(boxplot.stats(XNE_brain$OT)$out)), ], "OT")))


boxplot(XNE_brain$CB)
boxplot.stats(XNE_brain$CB)$out
outliers_data<-rbind(outliers_data, unlist(c(XNE_brain[which(XNE_brain$CB %in% c(boxplot.stats(XNE_brain$CB)$out)), ], "CB")))


boxplot(XNE_brain$HP)
boxplot.stats(XNE_brain$HP)$out
outliers_data<-rbind(outliers_data, unlist(c(XNE_brain[which(XNE_brain$HP %in% c(boxplot.stats(XNE_brain$HP)$out)), ], "HP")))
outliers_data<-rbind(outliers_data, unlist(c(XNE_brain[which(XNE_brain$HP %in% c(boxplot.stats(XNE_brain$HP)$out[2])), ], "HP")))

boxplot(XNE_brain$OB)
boxplot.stats(XNE_brain$OB)$out
outliers_data<-rbind(outliers_data, unlist(c(XNE_brain[which(XNE_brain$OB %in% c(boxplot.stats(XNE_brain$OB)$out)), ], "OB")))

boxplot(XNE_brain$DM)
boxplot.stats(XNE_brain$DM)$out
outliers_data<-rbind(outliers_data, unlist(c(XNE_brain[which(XNE_brain$DM %in% c(boxplot.stats(XNE_brain$DM)$out)), ], "DM")))


#XHE
XHE_brain<-brain_dataset[brain_dataset$Species=="XHE",]

boxplot(XHE_brain$TL)
boxplot.stats(XHE_brain$TL)$out

boxplot(XHE_brain$OT)
boxplot.stats(XHE_brain$OT)$out

boxplot(XHE_brain$CB)
boxplot.stats(XHE_brain$CB)$out

boxplot(XHE_brain$HP)
boxplot.stats(XHE_brain$HP)$out
outliers_data<-rbind(outliers_data, unlist(c(XHE_brain[which(XHE_brain$HP %in% c(boxplot.stats(XHE_brain$HP)$out[1])), ], "HP")))
outliers_data<-rbind(outliers_data, unlist(c(XHE_brain[which(XHE_brain$HP %in% c(boxplot.stats(XHE_brain$HP)$out[2])), ], "HP")))

boxplot(XHE_brain$OB)
boxplot.stats(XHE_brain$OB)$out
outliers_data<-rbind(outliers_data, unlist(c(XHE_brain[which(XHE_brain$OB %in% c(boxplot.stats(XHE_brain$OB)$out)), ], "OB")))

boxplot(XHE_brain$DM)
boxplot.stats(XHE_brain$DM)$out

## GVI
GVI_brain<-brain_dataset[brain_dataset$Species=="GVI",]

boxplot(GVI_brain$TL)
boxplot.stats(GVI_brain$TL)$out

boxplot(GVI_brain$OT)
boxplot.stats(GVI_brain$OT)$out

boxplot(GVI_brain$CB)
boxplot.stats(GVI_brain$CB)$out

boxplot(GVI_brain$HP)
boxplot.stats(GVI_brain$HP)$out
outliers_data<-rbind(outliers_data, unlist(c(GVI_brain[which(GVI_brain$HP %in% c(boxplot.stats(GVI_brain$HP)$out)), ], "HP")))

boxplot(GVI_brain$OB)
boxplot.stats(GVI_brain$OB)$out
outliers_data<-rbind(outliers_data, unlist(c(GVI_brain[which(GVI_brain$OB %in% c(boxplot.stats(GVI_brain$OB)$out)), ], "OB")))


boxplot(GVI_brain$DM)
boxplot.stats(GVI_brain$DM)$out


## XVA
XVA_brain<-brain_dataset[brain_dataset$Species=="XVA",]

boxplot(XVA_brain$TL)
boxplot.stats(XVA_brain$TL)$out
outliers_data<-rbind(outliers_data, unlist(c(XVA_brain[which(XVA_brain$TL %in% c(boxplot.stats(XVA_brain$TL)$out)), ], "TL")))
outliers_data<-rbind(outliers_data, unlist(c(XVA_brain[which(XVA_brain$TL %in% c(boxplot.stats(XVA_brain$TL)$out[2])), ], "TL")))

boxplot(XVA_brain$OT)
boxplot.stats(XVA_brain$OT)$out
outliers_data<-rbind(outliers_data, unlist(c(XVA_brain[which(XVA_brain$OT %in% c(boxplot.stats(XVA_brain$OT)$out)), ], "OT")))
outliers_data<-rbind(outliers_data, unlist(c(XVA_brain[which(XVA_brain$OT %in% c(boxplot.stats(XVA_brain$OT)$out[2])), ], "OT")))

boxplot(XVA_brain$CB)
boxplot.stats(XVA_brain$CB)$out
outliers_data<-rbind(outliers_data, unlist(c(XVA_brain[which(XVA_brain$CB %in% c(boxplot.stats(XVA_brain$CB)$out)), ], "CB")))

boxplot(XVA_brain$HP)
boxplot.stats(XVA_brain$HP)$out
outliers_data<-rbind(outliers_data, unlist(c(XVA_brain[which(XVA_brain$HP %in% c(boxplot.stats(XVA_brain$HP)$out)), ], "HP")))

boxplot(XVA_brain$OB)
boxplot.stats(XVA_brain$OB)$out
outliers_data<-rbind(outliers_data, unlist(c(XVA_brain[which(XVA_brain$OB %in% c(boxplot.stats(XVA_brain$OB)$out)), ], "OB")))


boxplot(XVA_brain$DM)
boxplot.stats(XVA_brain$DM)$out
outliers_data<-rbind(outliers_data, unlist(c(XVA_brain[which(XVA_brain$DM %in% c(boxplot.stats(XVA_brain$DM)$out)), ], "DM")))


## PQA
PQA_brain<-brain_dataset[brain_dataset$Species=="PQA",]

boxplot(PQA_brain$TL)
boxplot.stats(PQA_brain$TL)$out
outliers_data<-rbind(outliers_data, unlist(c(PQA_brain[which(PQA_brain$TL %in% c(boxplot.stats(PQA_brain$TL)$out)), ], "TL")))

boxplot(PQA_brain$OT)
boxplot.stats(PQA_brain$OT)$out
outliers_data<-rbind(outliers_data, unlist(c(PQA_brain[which(PQA_brain$OT %in% c(boxplot.stats(PQA_brain$OT)$out)), ], "OT")))

boxplot(PQA_brain$CB)
boxplot.stats(PQA_brain$CB)$out
outliers_data<-rbind(outliers_data, unlist(c(PQA_brain[which(PQA_brain$CB %in% c(boxplot.stats(PQA_brain$CB)$out)), ], "CB")))

boxplot(PQA_brain$HP)
boxplot.stats(PQA_brain$HP)$out
outliers_data<-rbind(outliers_data, unlist(c(PQA_brain[which(PQA_brain$HP %in% c(boxplot.stats(PQA_brain$HP)$out)), ], "HP")))

boxplot(PQA_brain$OB)
boxplot.stats(PQA_brain$OB)$out
outliers_data<-rbind(outliers_data, unlist(c(PQA_brain[which(PQA_brain$OB %in% c(boxplot.stats(PQA_brain$OB)$out)), ], "OB")))

boxplot(PQA_brain$DM)
boxplot.stats(PQA_brain$DM)$out


## NEL
NEL_brain<-brain_dataset[brain_dataset$Species=="NEL",]

boxplot(NEL_brain$TL)
boxplot.stats(NEL_brain$TL)$out
outliers_data<-rbind(outliers_data, unlist(c(NEL_brain[which(NEL_brain$TL %in% c(boxplot.stats(NEL_brain$TL)$out)), ], "TL")))

boxplot(NEL_brain$OT)
boxplot.stats(NEL_brain$OT)$out
outliers_data<-rbind(outliers_data, unlist(c(NEL_brain[which(NEL_brain$OT %in% c(boxplot.stats(NEL_brain$OT)$out)), ], "OT")))

boxplot(NEL_brain$CB)
boxplot.stats(NEL_brain$CB)$out
outliers_data<-rbind(outliers_data, unlist(c(NEL_brain[which(NEL_brain$CB %in% c(boxplot.stats(NEL_brain$CB)$out)), ], "CB")))

boxplot(NEL_brain$HP)
boxplot.stats(NEL_brain$HP)$out
outliers_data<-rbind(outliers_data, unlist(c(NEL_brain[which(NEL_brain$HP %in% c(boxplot.stats(NEL_brain$HP)$out)), ], "HP")))


boxplot(NEL_brain$OB)
boxplot.stats(NEL_brain$OB)$out

boxplot(NEL_brain$DM)
boxplot.stats(NEL_brain$DM)$out


## PMX
PMX_brain<-brain_dataset[brain_dataset$Species=="PMX",]

boxplot(PMX_brain$TL)
boxplot.stats(PMX_brain$TL)$out

boxplot(PMX_brain$OT)
boxplot.stats(PMX_brain$OT)$out

boxplot(PMX_brain$CB)
boxplot.stats(PMX_brain$CB)$out
outliers_data<-rbind(outliers_data, unlist(c(PMX_brain[which(PMX_brain$CB %in% c(boxplot.stats(PMX_brain$CB)$out)), ], "CB")))

boxplot(PMX_brain$HP)
boxplot.stats(PMX_brain$HP)$out

boxplot(PMX_brain$OB)
boxplot.stats(PMX_brain$OB)$out

boxplot(PMX_brain$DM)
boxplot.stats(PMX_brain$DM)$out
outliers_data<-rbind(outliers_data, unlist(c(PMX_brain[which(PMX_brain$DM %in% c(boxplot.stats(PMX_brain$DM)$out)), ], "DM")))


## LVI
LVI_brain<-brain_dataset[brain_dataset$Species=="LVI",]

boxplot(LVI_brain$TL)
boxplot.stats(LVI_brain$TL)$out
outliers_data<-rbind(outliers_data, unlist(c(LVI_brain[which(LVI_brain$TL %in% c(boxplot.stats(LVI_brain$TL)$out)), ], "TL")))

boxplot(LVI_brain$OT)
boxplot.stats(LVI_brain$OT)$out
outliers_data<-rbind(outliers_data, unlist(c(LVI_brain[which(LVI_brain$OT %in% c(boxplot.stats(LVI_brain$OT)$out)), ], "OT")))

boxplot(LVI_brain$CB)
boxplot.stats(LVI_brain$CB)$out
outliers_data<-rbind(outliers_data, unlist(c(LVI_brain[which(LVI_brain$CB %in% c(boxplot.stats(LVI_brain$CB)$out)), ], "CB")))

boxplot(LVI_brain$HP)
boxplot.stats(LVI_brain$HP)$out
outliers_data<-rbind(outliers_data, unlist(c(LVI_brain[which(LVI_brain$HP %in% c(boxplot.stats(LVI_brain$HP)$out)), ], "HP")))

boxplot(LVI_brain$OB)
boxplot.stats(LVI_brain$OB)$out
outliers_data<-rbind(outliers_data, unlist(c(LVI_brain[which(LVI_brain$OB %in% c(boxplot.stats(LVI_brain$OB)$out)), ], "OB")))

boxplot(LVI_brain$DM)
boxplot.stats(LVI_brain$DM)$out
outliers_data<-rbind(outliers_data, unlist(c(LVI_brain[which(LVI_brain$DM %in% c(boxplot.stats(LVI_brain$DM)$out[1])), ], "DM")))
outliers_data<-rbind(outliers_data, unlist(c(LVI_brain[which(LVI_brain$DM %in% c(boxplot.stats(LVI_brain$DM)$out[2])), ], "DM")))


## POB
POB_brain<-brain_dataset[brain_dataset$Species=="POB",]

boxplot(POB_brain$TL)
boxplot.stats(POB_brain$TL)$out

boxplot(POB_brain$OT)
boxplot.stats(POB_brain$OT)$out
outliers_data<-rbind(outliers_data, unlist(c(POB_brain[which(POB_brain$OT %in% c(boxplot.stats(POB_brain$OT)$out)), ], "OT")))

boxplot(POB_brain$CB)
boxplot.stats(POB_brain$CB)$out
outliers_data<-rbind(outliers_data, unlist(c(POB_brain[which(POB_brain$CB %in% c(boxplot.stats(POB_brain$CB)$out)), ], "CB")))

boxplot(POB_brain$HP)
boxplot.stats(POB_brain$HP)$out

boxplot(POB_brain$OB)
boxplot.stats(POB_brain$OB)$out

boxplot(POB_brain$DM)
boxplot.stats(POB_brain$DM)$out


## PGR
PGR_brain<-brain_dataset[brain_dataset$Species=="PGR",]

boxplot(PGR_brain$TL)
boxplot.stats(PGR_brain$TL)$out
outliers_data<-rbind(outliers_data, unlist(c(PGR_brain[which(PGR_brain$TL %in% c(boxplot.stats(PGR_brain$TL)$out)), ], "TL")))

boxplot(PGR_brain$OT)
boxplot.stats(PGR_brain$OT)$out
outliers_data<-rbind(outliers_data, unlist(c(PGR_brain[which(PGR_brain$OT %in% c(boxplot.stats(PGR_brain$OT)$out)), ], "OT")))


boxplot(PGR_brain$CB)
boxplot.stats(PGR_brain$CB)$out
outliers_data<-rbind(outliers_data, unlist(c(PGR_brain[which(PGR_brain$CB %in% c(boxplot.stats(PGR_brain$CB)$out)), ], "CB")))

boxplot(PGR_brain$HP)
boxplot.stats(PGR_brain$HP)$out

boxplot(PGR_brain$OB)
boxplot.stats(PGR_brain$OB)$out

boxplot(PGR_brain$DM)
boxplot.stats(PGR_brain$DM)$out


## PLA
PLA_brain<-brain_dataset[brain_dataset$Species=="PLA",]

boxplot(PLA_brain$TL)
boxplot.stats(PLA_brain$TL)$out

boxplot(PLA_brain$OT)
boxplot.stats(PLA_brain$OT)$out
outliers_data<-rbind(outliers_data, unlist(c(PLA_brain[which(PLA_brain$OT %in% c(boxplot.stats(PLA_brain$OT)$out)), ], "OT")))

boxplot(PLA_brain$CB)
boxplot.stats(PLA_brain$CB)$out
outliers_data<-rbind(outliers_data, unlist(c(PLA_brain[which(PLA_brain$CB %in% c(boxplot.stats(PLA_brain$CB)$out)), ], "CB")))

boxplot(PLA_brain$HP)
boxplot.stats(PLA_brain$HP)$out
outliers_data<-rbind(outliers_data, unlist(c(PLA_brain[which(PLA_brain$HP %in% c(boxplot.stats(PLA_brain$HP)$out)), ], "HP")))


boxplot(PLA_brain$OB)
boxplot.stats(PLA_brain$OB)$out
outliers_data<-rbind(outliers_data, unlist(c(PLA_brain[which(PLA_brain$OB %in% c(boxplot.stats(PLA_brain$OB)$out)), ], "OB")))


boxplot(PLA_brain$DM)
boxplot.stats(PLA_brain$DM)$out


## LPE
LPE_brain<-brain_dataset[brain_dataset$Species=="LPE",]

boxplot(LPE_brain$TL)
boxplot.stats(LPE_brain$TL)$out
outliers_data<-rbind(outliers_data, unlist(c(LPE_brain[which(LPE_brain$TL %in% c(boxplot.stats(LPE_brain$TL)$out)), ], "TL")))


boxplot(LPE_brain$OT)
boxplot.stats(LPE_brain$OT)$out
outliers_data<-rbind(outliers_data, unlist(c(LPE_brain[which(LPE_brain$OT %in% c(boxplot.stats(LPE_brain$OT)$out)), ], "OT")))


boxplot(LPE_brain$CB)
boxplot.stats(LPE_brain$CB)$out
outliers_data<-rbind(outliers_data, unlist(c(LPE_brain[which(LPE_brain$CB %in% c(boxplot.stats(LPE_brain$CB)$out)), ], "CB")))
# RUN only this part:  unlist(c(LPE_brain[which(LPE_brain$CB %in% c(boxplot.stats(LPE_brain$CB)$out)), ], "CB")) part to see what diss ID the outlier has

boxplot(LPE_brain$HP)
boxplot.stats(LPE_brain$HP)$out
outliers_data<-rbind(outliers_data, unlist(c(LPE_brain[which(LPE_brain$HP %in% c(boxplot.stats(LPE_brain$HP)$out)), ], "HP")))

boxplot(LPE_brain$OB)
boxplot.stats(LPE_brain$OB)$out
outliers_data<-rbind(outliers_data, unlist(c(LPE_brain[which(LPE_brain$OB %in% c(boxplot.stats(LPE_brain$OB)$out)), ], "OB")))

boxplot(LPE_brain$DM)
boxplot.stats(LPE_brain$DM)$out
outliers_data<-rbind(outliers_data, unlist(c(LPE_brain[which(LPE_brain$DM %in% c(boxplot.stats(LPE_brain$DM)$out)), ], "DM")))


## PWI
PWI_brain<-brain_dataset[brain_dataset$Species=="PWI",]

boxplot(PWI_brain$TL)
boxplot.stats(PWI_brain$TL)$out

boxplot(PWI_brain$OT)
boxplot.stats(PWI_brain$OT)$out

boxplot(PWI_brain$CB)
boxplot.stats(PWI_brain$CB)$out
outliers_data<-rbind(outliers_data, unlist(c(PWI_brain[which(PWI_brain$CB %in% c(boxplot.stats(PWI_brain$CB)$out)), ], "CB")))

boxplot(PWI_brain$HP)
boxplot.stats(PWI_brain$HP)$out

boxplot(PWI_brain$OB)
boxplot.stats(PWI_brain$OB)$out
outliers_data<-rbind(outliers_data, unlist(c(PWI_brain[which(PWI_brain$OB %in% c(boxplot.stats(PWI_brain$OB)$out)), ], "OB")))

boxplot(PWI_brain$DM)
boxplot.stats(PWI_brain$DM)$out
outliers_data<-rbind(outliers_data, unlist(c(PWI_brain[which(PWI_brain$DM %in% c(boxplot.stats(PWI_brain$DM)$out)), ], "DM")))


## PGR
PGR_brain<-brain_dataset[brain_dataset$Species=="PGR",]

boxplot(PGR_brain$TL)
boxplot.stats(PGR_brain$TL)$out
outliers_data<-rbind(outliers_data, unlist(c(PGR_brain[which(PGR_brain$TL %in% c(boxplot.stats(PGR_brain$TL)$out)), ], "TL")))

boxplot(PGR_brain$OT)
boxplot.stats(PGR_brain$OT)$out
outliers_data<-rbind(outliers_data, unlist(c(PGR_brain[which(PGR_brain$OT %in% c(boxplot.stats(PGR_brain$OT)$out)), ], "OT")))

boxplot(PGR_brain$CB)
boxplot.stats(PGR_brain$CB)$out
outliers_data<-rbind(outliers_data, unlist(c(PGR_brain[which(PGR_brain$CB %in% c(boxplot.stats(PGR_brain$CB)$out)), ], "CB")))

boxplot(PGR_brain$HP)
boxplot.stats(PGR_brain$HP)$out
outliers_data<-rbind(outliers_data, unlist(c(PGR_brain[which(PGR_brain$HP %in% c(boxplot.stats(PGR_brain$HP)$out)), ], "HP")))

boxplot(PGR_brain$OB)
boxplot.stats(PGR_brain$OB)$out

boxplot(PGR_brain$DM)
boxplot.stats(PGR_brain$DM)$out


## GME
GME_brain<-brain_dataset[brain_dataset$Species=="GME",]

boxplot(GME_brain$TL)
boxplot.stats(GME_brain$TL)$out
outliers_data<-rbind(outliers_data, unlist(c(GME_brain[which(GME_brain$TL %in% c(boxplot.stats(GME_brain$TL)$out)), ], "TL")))

boxplot(GME_brain$OT)
boxplot.stats(GME_brain$OT)$out
outliers_data<-rbind(outliers_data, unlist(c(GME_brain[which(GME_brain$OT %in% c(boxplot.stats(GME_brain$OT)$out)), ], "OT")))

boxplot(GME_brain$CB)
boxplot.stats(GME_brain$CB)$out
outliers_data<-rbind(outliers_data, unlist(c(GME_brain[which(GME_brain$CB %in% c(boxplot.stats(GME_brain$CB)$out)), ], "CB")))

boxplot(GME_brain$HP)
boxplot.stats(GME_brain$HP)$out
outliers_data<-rbind(outliers_data, unlist(c(GME_brain[which(GME_brain$HP %in% c(boxplot.stats(GME_brain$HP)$out)), ], "HP")))

boxplot(GME_brain$OB)
boxplot.stats(GME_brain$OB)$out
outliers_data<-rbind(outliers_data, unlist(c(GME_brain[which(GME_brain$OB %in% c(boxplot.stats(GME_brain$OB)$out)), ], "OB")))

boxplot(GME_brain$DM)
boxplot.stats(GME_brain$DM)$out
outliers_data<-rbind(outliers_data, unlist(c(GME_brain[which(GME_brain$DM %in% c(boxplot.stats(GME_brain$DM)$out)), ], "DM")))


## GAF
GAF_brain<-brain_dataset[brain_dataset$Species=="GAF",]

boxplot(GAF_brain$TL)
boxplot.stats(GAF_brain$TL)$out

boxplot(GAF_brain$OT)
boxplot.stats(GAF_brain$OT)$out

boxplot(GAF_brain$CB)
boxplot.stats(GAF_brain$CB)$out

boxplot(GAF_brain$HP)
boxplot.stats(GAF_brain$HP)$out

boxplot(GAF_brain$OB)
boxplot.stats(GAF_brain$OB)$out

boxplot(GAF_brain$DM)
boxplot.stats(GAF_brain$DM)$out


#create a table with outlier data

write.table(outliers_data,"outlier_data.txt", sep="\t", row.names = FALSE)



