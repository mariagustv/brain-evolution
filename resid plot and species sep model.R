
### Script for residual plots and for mcmc model separately for each species
## Author: Maria Gustavsson
## May 2025, (R version 4.4.3)

#Import packages
library(MCMCglmm)
library(dplyr)
library(ggplot2)


#Import the data from previous scripts 
brain_data <- read.table("brain_dataset.txt", header = T)
brain_data_sex <- read.table("brain_dataset_means_sex.txt", header = T) 


############################### Residual plots for females and males #############3

#add the abbriviations for the species to the new dataset
spp_match <- read.csv(file = "abb_table.csv", header=F, sep = ";")  #import table with species names list and abbreviations
#create new dataset called bd to save the abbriviations

bd <- brain_data_sex %>%
  left_join(spp_match, by = c("spp_phylo" = "V2"))

#### Plotting residuals for F and M with all species


# For Females 
female_data <- bd %>% filter(Sex == "F")
ggplot(female_data, aes(x = V1, y = brain_tot_resid)) +
  geom_bar(stat = "identity", fill = "red") +
  labs(title = "Total Brain residuals by Species (Females)", x = "Species")

# For Males
male_data <- bd %>% filter(Sex == "M")
ggplot(male_data, aes(x = V1, y = brain_tot_resid)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(title = "Total Brain residuals by Species (Males)", x = "Species")



################################################################################

### MCMC model for all of species run separately for each brain region


df <- brain_data #rename the dataframe to df. brain_data is the individual data (not log), and no means

species_list <- unique(df$Species) #Create a list with all unique species 
brain_regions <- c("TL", "OT", "CB", "HP", "OB", "DM", "brain_tot") #create a vector with the different brain regions
results_list <- list() #empty list to store the results



#Loop over the species
for (sp in species_list) {
  df_sp <- df %>% filter(Species == sp)
  #loop over the brain regions  
  for (brain in brain_regions) {
    # Only fit the model if there is enough data and no missing values 
    df_sub <- df_sp %>% filter(!is.na(.data[[brain]]), !is.na(Weight))
    
    if (nrow(df_sub) > 2) { #First this was run with >5 to inlcude only the species with over 5 indiviudals. Reduced later to 2 to inlcude more species but they must be interpreted with caution since they have small sample size
      prior <- list(R = list(V = 1, nu = 0.002)) #weak informed priors
      
      formula <- as.formula(paste(brain, "~Weight")) #Create the model 
      
      model <- tryCatch(
        MCMCglmm(
          formula,
          data = df_sub,
          family = "gaussian",
          prior = prior,
          nitt = 13000,
          burnin = 3000,
          thin = 10
        ),
        error = function(e) {
          print(paste("Model error for", sp, brain, ":", e$message)) #print error message if there is an error
          NULL
        }
      )
      #extract and store the results in the list
      if (!is.null(model)) {
        summary_stats <- summary(model)$solutions
        results_list[[paste(sp, brain, sep = "_")]] <- data.frame(
          Species = sp,
          brain_regions = brain,
          mean_estimate = summary_stats["Weight", "post.mean"],
          l95 = summary_stats["Weight", "l-95% CI"],
          u95 = summary_stats["Weight", "u-95% CI"],
          pMCMC = summary_stats["Weight", "pMCMC"]
        )
      }
    }
  }
}
# all results are combined into one variable called final_results
final_results <- do.call(rbind, results_list[!sapply(results_list, is.null)])
print(final_results)

## some species has only one or two individuals and is therefore not included in the model (HMI, GSE, PPR, LVE)
#The species that has less than/or 5 individuals and might be reliable: LME, LNI, HFO, LDO, GAF, XMY
##################################################################################

## Count the number of females and males for each species
brain_data %>%
  group_by(Species, Sex) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(Species, desc(count))

sum_table <- brain_data %>%
  group_by(Species, Sex, Origin) %>%
  summarise(count = n(), .groups = "drop")
  arrange(Species, Sex, Origin)

print(sum_table)




