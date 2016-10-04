# GK Fragiadakis
# December 8th, 2015
# Processing neonatal data (QBaca)

library(dplyr)
library(blackbuck)
library(FDAlibrary)
library(ggplot2)
library(samr)

#-------------Signaling Data Structures----------------#

dataDirectory <- "/Users/gabrielafragiadakis/Documents/neonatal/exported_statistics/"
raw_df <- read.csv(paste(dataDirectory, "neonatals-medians.csv", sep = ""))

raw_df <- raw_df %>%
  dplyr::select(-(FCS.Filename)) %>%
  dplyr::rename(Pair = Individuals) %>%
  dplyr::mutate(Individual = as.factor(paste(Pair, SampleType, sep = "_"))) %>%
  dplyr::rename(Condition = Conditions)

# Get folds
fold_df <- get_folds(raw_df, basal_name = "unstim", fold_type = "asinh", unique_IDs = "filename",
                     annotation_columns = c("Individual", "Pair","SampleType", "Condition"))

# Unstim data
unstim_df <- dplyr::filter(raw_df, Condition == "unstim")

#transform unstim data
unstim_df <- reshape2::melt(unstim_df)
unstim_df$value <- asinh(unstim_df$value / 5)
unstim_df <- reshape2::dcast(unstim_df, ... ~ variable)
unstim_df <- dplyr::select(unstim_df, -filename)
Unstim_df <- unstim_df


# building blocks: unstim_df, LPS_df, Cocktail_df
# Add condition names to features
annotation_columns = c("Individual", "Pair","SampleType", "Condition")

LPS_df <- dplyr::filter(fold_df, Condition == "LPS")
colnames(LPS_df) <- paste(colnames(LPS_df), "_LPS", sep = "")
colnames(LPS_df)[1:4] <- annotation_columns
Cocktail_df <- dplyr::filter(fold_df, Condition == "Cocktail")
colnames(Cocktail_df) <- paste(colnames(Cocktail_df), "_Cocktail", sep = "")
colnames(Cocktail_df)[1:4] <- annotation_columns
new_df <- cbind(LPS_df, Cocktail_df[,!(names(Cocktail_df) %in% annotation_columns)])
  


##-----------abundance data----------## 

abundance_df <- read.csv(paste(dataDirectory, "neonatals-abundances.csv", sep = ""))

abundance_df <- abundance_df %>%
  dplyr::select(-(FCS.Filename)) %>%
  dplyr::rename(Pair = Individuals) %>%
  dplyr::mutate(Individual = as.factor(paste(Pair, SampleType, sep = "_"))) %>%
  dplyr::rename(Condition = Conditions)  %>%
  dplyr::filter(Condition == "unstim")

base_counts <- abundance_df$Mononuclearcells_EventCounts
singlet_counts <- abundance_df$Singlets_EventCounts
cell_types <- abundance_df[,!(names(abundance_df) %in% c(annotation_columns, 
                                                         "Singlets_EventCounts", "Mononuclearcells_EventCounts", 
                                                         "Granulocytes_EventCounts", "filename"))]
grans <- abundance_df$Granulocytes_EventCounts
percents <- cbind(cell_types/base_counts, grans/singlet_counts)
percents_df <- data.frame(filename = abundance_df$filename, abundance_df[, colnames(abundance_df) %in% annotation_columns], percents)
# write.csv(percents_df, file = "~/Documents/neonatal/neonatals_percentages.csv")


