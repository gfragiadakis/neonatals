library(devtools)
install_github("nolanlab/scaffold", ref = "multiFilesClustering")
library(scaffold)
library(FDAlibrary)
library(dplyr)
library(blackbuck)
library(samr)
library(reshape2)
library(ggplot2)

annotation_columns = c("Individual", "Pair","SampleType", "Condition")


cluster_full <- read.table("~/Documents/neonatal/scaffold/200_clusters_new/scaffold1_file/030715_UC_Neonatal_UC14_M_Unstim_M_unstim_Singlets.fcs.clustered.txt",
                           sep = "\t", check.names = FALSE)
proteins <- c("CREB","MAPKAPK2","NFkB", "p38", "IkB", "S6","STAT1","STAT3","STAT5", "ERK1-2")

# 200 x 54 table of cluster medians overall (here is where I will add the coloring index)
cluster_medians <- cluster_full[2:201, 1:54]
colnames(cluster_medians) <- sapply(cluster_full[1, 1:54], as.character)

# cluster_stats has all the info from the individual files
cluster_stats <- cluster_full[2:201, -c(1:54)]
colnames(cluster_stats) <- sapply(cluster_full[1, -c(1:54)], as.character)
rownames(cluster_stats) <- paste("@cluster", 1:200, sep = "_")

# get file list
files <- colnames(cluster_stats)
splitted <- strsplit(files, split = "@")
files <- lapply(splitted, "[[", 2)
file_list <- unique(files)

# loop through
data_matrix <- mat.or.vec(nr = 60, nc = 2000)
rownames(data_matrix) <- unlist(file_list)

for (i in 1:length(file_list)){
  file <- file_list[[i]]
  ind <- grep(pattern = file, colnames(cluster_stats))
  feat_vals <- cluster_stats[, ind]
  new_names <- strsplit(colnames(feat_vals), split = "@")
  new_names <- sapply(new_names, "[[", 1)
  colnames(feat_vals) <- new_names
  values <- feat_vals[, colnames(feat_vals) %in% proteins]
  colnames(values) <- paste(colnames(values), "_", sep = "")
  value_vec <- unlist(values)
  # print(names(value_vec) == colnames(data_matrix))
  colnames(data_matrix) <- names(value_vec)
  value_vec_numeric <- as.numeric(as.character(value_vec))
  data_matrix[i, ] <- value_vec_numeric
  # print(file == rownames(data_matrix)[i])
}

# Ok have data matrix with everything, like raw_df, now just need annotation
filenames <- rownames(data_matrix)
split_files <- strsplit(filenames, split = "_")
SampleType <- sapply(split_files, "[[", 5)
Condition <- sapply(split_files, "[[", 6)
Pair <- sapply(split_files, "[[", 4)
Individual <- paste(Pair, SampleType, sep = "_")

cluster_df <- data.frame(filename = filenames, 
                         Condition = Condition, 
                         SampleType = SampleType,
                         Pair = Pair, 
                         Individual = Individual, 
                         data_matrix)

# get ones per condition: 
fold_df <- get_folds(cluster_df, basal_name = "Unstim", fold_type = "asinh", unique_IDs = "filename",
                     annotation_columns = annotation_columns)
LPS_df <- dplyr::filter(fold_df, Condition == "LPS")
Cocktail_df <- dplyr::filter(fold_df, Condition == "Cocktail")


unstim_df <- dplyr::filter(cluster_df, Condition == "Unstim")
unstim_df <- reshape2::melt(unstim_df)
unstim_df$value <- asinh(unstim_df$value / 5)
unstim_df <- reshape2::dcast(unstim_df, ... ~ variable)
unstim_df <- dplyr::select(unstim_df, -filename)
Unstim_df <- unstim_df


# now have Unstim_df, Cocktail_df, LPS_df

y <- dplyr::select(Unstim_df, Pair, SampleType)
y <- cbind(y, resp = rep(1:10, each  = 2))
for (i in 1:nrow(y)){
  if ( y[i, "SampleType"] == "B"){
    y[i, "resp"] <- -(y[i, "resp"])
  }
}
y <- y$resp

dfs <- list("LPS" = LPS_df, "Cocktail" = Cocktail_df, "Unstim" = Unstim_df)
neonatal_SAM <- function(dfs, y, resp.type, SAMDirectory){
  output <- list()
  for (i in 1:length(dfs)){
    df <- dfs[[i]]
    x <- df[,!(names(df) %in% annotation_columns)]
    sink("/dev/null")
    model <-  SAM(t(x), y, genenames = colnames(x), resp.type = resp.type, nperms = 1000, fdr.output = 0.01)
    sink()
    sig.up <- model$siggenes.table$genes.up[model$siggenes.table$genes.up[,"q-value(%)"] == "0", ]
    
    
    sig.lo <- model$siggenes.table$genes.lo[model$siggenes.table$genes.lo[,"q-value(%)"] == "0", ]
    
    sigs <- rbind(sig.up, sig.lo)
    # write.csv(sigs, file = paste(SAMDirectory, names(dfs)[i], "_fdr_0.01.csv", sep = ""))
    output[[names(dfs)[i]]] <- sigs
    
    sig_df <- data.frame(sigs)
    feat_vals <- cbind(Individual = df$Individual, df[, colnames(df) %in% sig_df$Gene.ID])
    # write.csv(feat_vals, file = paste(SAMDirectory, names(dfs)[i], "_feature_values_fdr_0.01.csv", sep = ""))
    
    df <- melt(df)
    df <- df[df[, "variable"] %in% sig_df$Gene.ID,]
    p <- ggplot(df, aes(variable, value)) + geom_boxplot(aes(colour = SampleType)) + ggtitle(names(dfs)[i]) + coord_flip() + theme(axis.text = element_text(size = 5))
    print(p)
    # ggsave(filename = paste(SAMDirectory, names(dfs)[i], "_sig_plot.pdf", sep = ""), width = 30, height = 20)
    print(names(dfs)[i])
    print("Higher in Mother")
    print(sig.up[, "Gene ID"])
    print("Higher in Baby")
    print(sig.lo[, "Gene ID"])
  }
  return(output)
}

paired_output <- neonatal_SAM(dfs, y, resp.type = "Two class paired", SAMDirectory = "./SAM/")

# here we want to find them, get the fold change, and filter
paired_output_thresh <- list()
sig_LPS <- dfs$LPS[, colnames(dfs$LPS) %in% paired_output$LPS[, "Gene ID"]]
thresh_LPS <- colnames(sig_LPS)[abs(apply(sig_LPS, 2, mean, na.rm = TRUE)) > 0.1]
paired_output_thresh[["LPS"]] <- paired_output$LPS[paired_output$LPS[, "Gene ID"] %in% thresh_LPS, ]

sig_Unstim <- dfs$Unstim[, colnames(dfs$Unstim) %in% paired_output$Unstim[, "Gene ID"]]
thresh_Unstim <- colnames(sig_Unstim)[abs(apply(sig_Unstim, 2, mean, na.rm = TRUE)) > 0.3]
paired_output_thresh[["Unstim"]] <- paired_output$Unstim[paired_output$Unstim[, "Gene ID"] %in% thresh_Unstim, ]

sig_Cocktail <- dfs$Cocktail[, colnames(dfs$Cocktail) %in% paired_output$Cocktail[, "Gene ID"]]
thresh_Cocktail <- colnames(sig_Cocktail)[abs(apply(sig_Cocktail, 2, mean, na.rm = TRUE)) > 0.3]
paired_output_thresh[["Cocktail"]] <- paired_output$Cocktail[paired_output$Cocktail[, "Gene ID"] %in% thresh_Cocktail, ]



# take paired_output, for each element, have a list of the Gene ID and the Score
full_dat <- c()
for (i in 1:3){
  # dat <- cbind(paired_output[[i]][, "Gene ID"], paired_output[[i]][, "Score(d)"])
  dat <- cbind(paired_output_thresh[[i]][, "Gene ID"], paired_output_thresh[[i]][, "Score(d)"])
  full_dat <- rbind(full_dat, dat)
}

split_hits <- strsplit(full_dat[, 1], split = "_")
protein_list <- sapply(split_hits, "[[", 1)
cluster_list <- sapply(split_hits, "[[", 2)
sig_data <- data.frame(protein = protein_list, cluster = cluster_list, direction = sign(as.numeric(full_dat[,2])))

SAM_count <- mat.or.vec(nr = 200, nc = 2)
colnames(SAM_count) <- c("cluster_number", "hits")
SAM_count[1:200, 1] <- c(1:200)

library(dplyr)
for (i in 1:200){
  # print(length(unique(dplyr::filter(sig_data,i cluster == i)[, 3])))
  total <- sum(dplyr::filter(sig_data, cluster == i)[, 3])
  SAM_count[i, 2] <- total
}


