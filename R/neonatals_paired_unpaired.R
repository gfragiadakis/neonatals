# # have the list of significant "genes" in the list output
# paired_output <- output
# 
# # We want to put in the original data frames and the SAM output, and get out the features from
# # those data frame and a list of features
# # start with output list
# 
# get_SAM_lists <- function(dfs, output){
#   feature_names <- lapply(output, function(x) x[, "Gene ID"])
#   
#   feature_values <- list()
#   for (i in names(dfs)){
#     df <- dfs[[i]]
#     feature_values[[i]] <- cbind(Individual = df[, "Individual"], df[, colnames(df) %in% feature_names[[i]]])
#   }
#   return(list(feature_names, feature_values))
# }
# 
# paired <- get_SAM_lists(dfs = dfs, output = paired_output)
# unpaired <- get_SAM_lists(dfs = dfs, output = unpaired_output)
# 
# paired[[2]] <- lapply(paired[[2]], function(x) dplyr::select(x, -Individual))
# 
# paired_x <- as.matrix(do.call(cbind, args = paired[[2]]))
# paired_df <- data.frame(Individual = dfs$Unstim$Individual, SampleType = dfs$Unstim$SampleType, Pair= dfs$Unstim$Pair, paired_x)
# 
# unpaired_x <- as.matrix(do.call(cbind, args = unpaired[[2]]))
# 
# ###--------PCA ------###
# # x <- paired_x
# x <- unpaired_x
# pr.out <- prcomp(x, scale = FALSE)
# transformed_x <- pr.out$x
# PC_df <- data.frame(transformed_x, Individual = paired_df$Individual, SampleType = paired_df$SampleType, Pair = paired_df$Pair)
# 
# 
# # loadings
# loadings1 <- pr.out$rotation[, "PC1"]
# loadings2 <- pr.out$rotation[, "PC2"]
# 
# # plotting
# pca <- ggplot(PC_df, aes(PC1, PC2))
# pca + geom_point(size = 3, aes(colour = SampleType)) + scale_color_manual(values=c("firebrick", "blue")) + ggtitle("PCA Unpaired")
# pca + geom_point(size = 3, aes(colour = Pair)) + ggtitle("PCA Unpaired")
# pca + geom_point(size = 3, aes(colour = Pair)) + ggtitle("PCA paired")


-----------------
# input df
# For signaling mom_df, baby_df
pairedDirectory <- "~/Documents/neonatal/paired/"
  
dfs <- list("LPS" = LPS_df, "Cocktail" = Cocktail_df, "Unstim" = Unstim_df)

df <- do.call(cbind, args = dfs)
mom_df <- df %>% dplyr::filter(Unstim.SampleType == "M") %>% dplyr::select(-contains("Individual"), -contains("Pair"), - contains("SampleType"), -contains("Condition"))
baby_df <- df %>% dplyr::filter(Unstim.SampleType == "B") %>%  dplyr::select(-contains("Individual"), -contains("Pair"), - contains("SampleType"), -contains("Condition"))

# for abundance version
# 
# mom_df <- percents_df %>% dplyr::filter(SampleType == "M") %>% dplyr::select(-filename, -Pair, -SampleType, -Condition, -Individual)
# baby_df <- percents_df %>% dplyr::filter(SampleType == "B") %>% dplyr::select(-filename, -Pair, -SampleType, -Condition, -Individual)

paired_pvals <- mapply(function(x, y) t.test(x,y, paired = TRUE)$p.value, mom_df, baby_df)
unpaired_pvals <- mapply(function(x, y) t.test(x,y, paired = FALSE)$p.value, mom_df, baby_df)

plot(-log(unpaired_pvals), -log(paired_pvals))

# Number of comparisons: 20 cell types x 10 phosphos x 3 conditions
bonf <- .05/(ncol(mom_df))

significance <- vector(mode = "character", length = length(unpaired_pvals))
significance[(unpaired_pvals < bonf) & (paired_pvals < bonf)] <- "Both"
significance[(unpaired_pvals < bonf) & (paired_pvals > bonf)] <- "Unpaired only"
significance[(unpaired_pvals > bonf) & (paired_pvals < bonf)] <- "Paired only"
significance[(unpaired_pvals > bonf) & (paired_pvals > bonf)] <- "Not significant"

plot_df <- data.frame(unpaired_pvals, paired_pvals, Significance = significance, Features = names(paired_pvals))
p <- ggplot(plot_df, aes(-log(unpaired_pvals), -log(paired_pvals))) + 
  geom_point(size = 1, aes(colour = Significance))  + 
  scale_color_manual(values=c("blue", " dark gray", "red", "orange"))
ggsave(paste(pairedDirectory, "scatter_paired_unpaired.pdf", sep = ""))

# we would take the ones that are significant, get the values, and spit that into the PCA

x <- df[, colnames(df) %in% names(paired_pvals[paired_pvals < bonf])]
x <- df[, colnames(df) %in% names(unpaired_pvals[unpaired_pvals < bonf])]

# WE JUST WANT TO USE THE ENTIRE X NOT A SUBSET

df_num <- df %>% dplyr::select(-contains("Individual"), -contains("Pair"), - contains("SampleType"), -contains("Condition"))
x <- as.matrix(df_num)

pr.out <- prcomp(x, scale = FALSE)
transformed_x <- pr.out$x
PC_df <- data.frame(transformed_x, Individual = unstim_df$Individual, SampleType = unstim_df$SampleType, Pair = unstim_df$Pair)
loadings1 <- pr.out$rotation[, "PC1"]
loadings2 <- pr.out$rotation[, "PC2"]

# plotting
pca <- ggplot(PC_df, aes(PC1, PC2, group = Pair))
pca + geom_point(size = 2, aes(colour = SampleType)) + scale_color_manual(values=c("firebrick", "blue")) + ggtitle("PCA by population")
ggsave(paste(pairedDirectory, "PCA_by_population.pdf"))

pca + geom_point(size = 2, aes(colour = Pair)) + ggtitle("PCA by mom-baby pair") + geom_line(aes(colour = Pair))
ggsave(paste(pairedDirectory, "PCA_by_pair.pdf"))

pca <- ggplot(PC_df, aes(Pair, PC2))
pca + geom_point(size = 2, aes(colour = Pair)) + ggtitle("PCA paired")
ggsave(paste(pairedDirectory, "PCA_by_pair_2.pdf"))


Unpaired_vals <- dplyr::filter(plot_df, Significance == "Unpaired only")
