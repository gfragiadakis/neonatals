# GK Fragiadakis
# December 8th, 2015
# Processing neonatal data (QBaca)

library(dplyr)
library(FDAlibrary)
library(ggplot2)
library(samr)

# read in signaling data
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


annotation_columns = c("Individual", "Pair","SampleType", "Condition")


# Add conditions onto it
LPS_df <- dplyr::filter(fold_df, Condition == "LPS")
colnames(LPS_df) <- paste(colnames(LPS_df), "_LPS", sep = "")
colnames(LPS_df)[1:4] <- annotation_columns
Cocktail_df <- dplyr::filter(fold_df, Condition == "Cocktail")
colnames(Cocktail_df) <- paste(colnames(Cocktail_df), "_Cocktail", sep = "")
colnames(Cocktail_df)[1:4] <- annotation_columns
new_df <- cbind(LPS_df, Cocktail_df[,!(names(Cocktail_df) %in% annotation_columns)])

# now have new_df for fold changed stims, and unstim_df for unstim.  Should run SAM two-class paired on both.  
# can also do unstim_df, LPS_df, Cocktail_df

# abundance data

abundance_df <- read.csv(paste(dataDirectory, "neonatals-abundances.csv", sep = ""))

abundance_df <- abundance_df %>%
  dplyr::select(-(FCS.Filename)) %>%
  dplyr::rename(Pair = Individuals) %>%
  dplyr::mutate(Individual = as.factor(paste(Pair, SampleType, sep = "_"))) %>%
  dplyr::rename(Condition = Conditions)  %>%
  dplyr::filter(Condition == "unstim")

base_counts <- abundance_df$Mononuclearcells_EventCounts
singlet_counts <- abundance_df$Singlets_EventCounts
cell_type_cols <- colnames(abundance_df)[7:23]
cell_types <- abundance_df[, colnames(abundance_df) %in% cell_type_cols]
grans <- abundance_df$Granulocytes_EventCounts
percents <- cbind(cell_types/base_counts, grans/singlet_counts)
percents_df <- data.frame(abundance_df[,c(1:4, 25)], percents)
write.csv(percents_df, file = "~/Documents/neonatal/neonatals_percentages.csv", rownames = FALSE)


