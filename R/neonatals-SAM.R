# GKFragiadakis
# Neonatal SAM analysis 
# data structures generated from neonatals-data-structures.R: unstim_df, fold_df

library(dplyr)
library(samr)
library(reshape2)

SAMDirectory <- "~/Documents/neonatal/SAM/"
# generate response varaible y for Two Class Paired (Mom and Baby)
y <- dplyr::select(unstim_df, Pair, SampleType)
y <- cbind(y, resp = rep(1:10, each  = 2))

for (i in 1:nrow(y)){
  if ( y[i, "SampleType"] == "B"){
    y[i, "resp"] <- -(y[i, "resp"])
  }
}

y <- y$resp
annotation_columns = c("Individual", "Pair","SampleType", "Condition")

LPS_df <- dplyr::filter(fold_df, Condition == "LPS")
Cocktail_df <- dplyr::filter(fold_df, Condition == "Cocktail")
Unstim_df <- unstim_df

dfs <- list("LPS" = LPS_df, "Cocktail" = Cocktail_df, "Unstim" = Unstim_df)

neonatal_SAM <- function(dfs, y, resp.type, SAMDirectory){
  output <- list()
  for (i in 1:length(dfs)){
    df <- dfs[[i]]
    x <- df[,!(names(df) %in% annotation_columns)]
    model <-  SAM(t(x), y, genenames = colnames(x), resp.type = resp.type, nperms = 1000, fdr.output = 0.01)
    sig.up <- model$siggenes.table$genes.up[model$siggenes.table$genes.up[,"q-value(%)"] == "0", ]
    sig.lo <- model$siggenes.table$genes.lo[model$siggenes.table$genes.lo[,"q-value(%)"] == "0", ]
    sigs <- rbind(sig.up, sig.lo)
    write.csv(sigs, file = paste(SAMDirectory, names(dfs)[i], "_fdr_0.01.csv", sep = ""))
    output[[names(dfs)[i]]] <- sigs
    
    sig_df <- data.frame(sigs)
    feat_vals <- cbind(Individual = df$Individual, df[, colnames(df) %in% sig_df$Gene.ID])
    write.csv(feat_vals, file = paste(SAMDirectory, names(dfs)[i], "_feature_values_fdr_0.01.csv", sep = ""))
    
    df <- melt(df)
    
    df <- df[df[, "variable"] %in% sig_df$Gene.ID,]
    ggplot(df, aes(variable, value)) + geom_boxplot(aes(colour = SampleType)) + theme(axis.text = element_text(size = rel(1), angle = 90)) + ggtitle(names(dfs)[i])
    ggsave(filename = paste(SAMDirectory, names(dfs)[i], "_sig_plot.pdf", sep = ""), width = 30, height = 20)
  }
  return(output)
}

paired_output <- neonatal_SAM(dfs, y, resp.type = "Two class paired", SAMDirectory <- "~/Documents/neonatal/SAM/")
# unpaired_output <- neonatal_SAM(dfs, y = ifelse(sign(y) > 0, 2, 1), resp.type = "Two class unpaired", SAMDirectory <- "~/Documents/neonatal/SAM/SAM_paired_v_unpaired/")


## Abundances! % mononuclear
# generated percents_df (includes annotations) and percents (numeric only) from neonatals data structures
# same y as above
x <- percents
model <-  SAM(t(x), y, genenames = colnames(x), resp.type = "Two class paired", nperms = 1000, fdr.output = 0.01)
sig.up <- model$siggenes.table$genes.up[model$siggenes.table$genes.up[,"q-value(%)"] == "0", ]
sig.lo <- model$siggenes.table$genes.lo[model$siggenes.table$genes.lo[,"q-value(%)"] == "0", ]
sigs <- rbind(sig.up, sig.lo)
write.csv(sigs, file = paste(SAMDirectory,  "percent_mononuclear_fdr_0.01.csv", sep = ""))

df <- percents_df
sig_df <- data.frame(sigs)
feat_vals <- cbind(Individual = df$Individual, df[, colnames(df) %in% sig_df$Gene.ID])
write.csv(feat_vals, file = paste(SAMDirectory, "percent_mononuclear_feature_values_fdr_0.01.csv", sep = ""))

df <- melt(df)

df <- df[df[, "variable"] %in% sig_df$Gene.ID,]
ggplot(df, aes(variable, value)) + geom_boxplot(aes(colour = SampleType)) + theme(axis.text = element_text(size = rel(1), angle = 90)) + ggtitle("Abundance differences (% Mononuclear cells)")
ggsave(filename = paste(SAMDirectory, "percent_mononuclear_sig_plot.pdf", sep = ""), width = 30, height = 20)

# plot all
df <- percents_df


df <- reshape2::melt(df)

ggplot(df, aes(variable, value)) + geom_boxplot(aes(colour = SampleType)) + theme(axis.text = element_text(size = rel(1), angle = 90)) + ggtitle("Abundance differences (% Mononuclear cells)")
ggsave(filename = paste(SAMDirectory, "percent_mononuclear_all_plot.pdf", sep = ""), width = 30, height = 20)

