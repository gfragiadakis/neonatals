## Script to correlate features in the mother to features in the baby
# December 15, 2015
# Author GKFragiadakis

#####-----Preprocessing (after running neonatals-data-structures.R)

# separate mom and baby
df <- new_df

separate_mvb <- function(df, annotation_columns, rename = TRUE){
  
  baby_df <- dplyr::filter(df, SampleType == "B")
  mother_df <- dplyr::filter(df, SampleType == "M")
  
  # test pairs line up
  if (identical(baby_df$Pair, mother_df$Pair)){
    
    # restrict to numeric columns
    baby_df <- baby_df[, !(colnames(baby_df) %in% annotation_columns)]
    mother_df <- mother_df[, !(colnames(mother_df) %in% annotation_columns)]
    
    if (rename == TRUE){
      colnames(baby_df) <- paste(colnames(baby_df), "B", sep = "_")
      colnames(mother_df) <- paste(colnames(mother_df), "M", sep = "_") 
    }
  df_list <- list(Baby = baby_df, Mother = mother_df)
  return(df_list)
  
  } else {
    print("Error: ordering doesn't match")
  }
  
}

output_directory <- "~/Documents/neonatal/correlation/"
dfs <- separate_mvb(new_df, annotation_columns, rename = TRUE)
baby_df <- dfs$Baby
mother_df <- dfs$Mother

####-------Correlation-------------
plot_correlations_asymmetric <- function(df1, df2, cor_threshold = 0.5, output_directory, background = "white", main_title = "basic"){
  
  # cluster rows and columns
  cormat <- cor(df1, df2)
  cormat_clustered <- cluster_matrix(cormat)
  clust <- cormat_clustered$data
  cormat_melted <- reshape2::melt(clust)
  
  # plot correlations
  cor_plot <- qplot(x=Var1, y=Var2, data=cormat_melted, fill=value, geom="tile") +
    scale_fill_gradient2(low = "red", high = "blue", mid = background, midpoint = 0, limits=c(-1, 1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
   ggsave(paste(output_directory, main_title, "correlation_map.pdf", sep = ""), plot = cor_plot, width = 40, height = 40)
  
  adj <- make_adjacency(clust, cor_threshold)
  adj_melted <- reshape2::melt(adj)
  
  cor_plot <- qplot(x=Var1, y=Var2, data=adj_melted, fill=value, geom="tile") +
    scale_fill_gradient2(low = "red", high = "blue", mid = background, midpoint = 0, limits=c(-1, 1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste(output_directory, main_title, "adj_map.pdf", sep = ""), plot = cor_plot, width = 40, height = 40)
}

plot_correlations_asymmetric(df1 = baby_df, df2 = mother_df, cor_threshold = 0.5, output_directory, background = "white", main_title = "Full_MvB_")

# conditions alone

df_conditions <- list(unstim_df, Cocktail_df, LPS_df)
df_names <- c("Unstim", "Cocktail", "LPS")

for (i in 1:length(df_names)){
  dfs <- separate_mvb(df_conditions[[i]], annotation_columns, rename = TRUE)
  df1 <- dfs$Baby
  df2 <- dfs$Mother
  plot_correlations_asymmetric(df1, df2, cor_threshold = 0.5, output_directory, background = "white", main_title = paste(df_names[i], "_MvB_",sep=""))
  
}


# R values
R_values <- sort(cormat_melted$value)
hist(R_values)
qplot(value, data = cormat_melted)
ggsave(paste(output_directory, "correlation_histogram.pdf", sep = ""))
write.csv(cormat_melted, paste(output_directory, "fold_data_correlations.csv", sep = ""))

######--------Adjacencies

cormat <- cor(baby_df, mother_df)
cormat_clustered <- cluster_matrix(cormat)
clust <- cormat_clustered$data
thresholds <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8)

for (i in thresholds){
  
  adj <- make_adjacency(clust, cor_threshold = i)
  adj_melted <- reshape2::melt(adj)
  
  cor_plot <- qplot(x=Var1, y=Var2, data=adj_melted, fill=value, geom="tile") +
              scale_fill_gradient2(low = "red", high = "blue", mid = "white", midpoint = 0, limits=c(-1, 1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  output_directory <- "~/Documents/neonatal/correlation/"
  ggsave(paste(output_directory, i, "adj_map_MvB.pdf", sep = ""), plot = cor_plot, width = 40, height = 40)
}


# diagonals

co_of_var <- function(x){
  return(sd(x)/mean(x))
}

baby_cov <- apply(baby_df, 2, co_of_var)
mother_cov <- apply(mother_df, 2, co_of_var)
covs <- cbind(baby_cov, mother_cov)
Rvalues <- diag(cormat)
diagonal_cors <- as.data.frame(cbind(covs, Rvalues))
sorted_cors_coeff <- diagonal_cors[order(diagonal_cors$Rvalues, decreasing = TRUE), ]


baby_var <- apply(baby_df, 2, var)
mother_var <- apply(mother_df, 2, var)
vars <- cbind(baby_var, mother_var)
Rvalues <- diag(cormat)
diagonal_cors <- data.frame(Baby_Feature = names(baby_var), Mother_Feature = names(mother_var), baby_var, mother_var, Rvalues)
sorted_cors_variance <- diagonal_cors[order(diagonal_cors$Rvalues, decreasing = TRUE), ]

write.csv(sorted_cors_coeff, file = paste(output_directory, "diagonals_cor_coefficient_of_vars.csv", sep = ""))
write.csv(sorted_cors_variance, file = paste(output_directory, "diagonals_cor_variance.csv", sep = ""))

# Values above a certain R

restricted <- diagonal_cors[diagonal_cors[, "Rvalues"] > 0.6, ]
restricted <- restricted[order(restricted$baby_var, decreasing = TRUE), ]

df_total <- cbind(baby_df, mother_df)
out_dir <- "~/Documents/neonatal/correlation/restricted_R_0.6/"

for (i in 1:nrow(restricted)){
  x <-as.character(restricted$Mother_Feature[i])
  y <- as.character(restricted$Baby_Feature[i])
  rplot <- ggplot(df_total, aes_string(x, y)) + geom_point()  + stat_smooth(method = "lm", se = FALSE)
  ggsave(paste(out_dir, i, "_R_", round(restricted$Rvalues[i], 2), "_plot.pdf", sep = ""), plot = rplot, width = 8, height =8)
}

