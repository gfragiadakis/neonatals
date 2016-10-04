naive_CD8 <- abundance_df$NaiveCD8.TcellsCD45RA._EventCounts/abundance_df$CD8.Tcell_EventCounts
memory_CD8 <- abundance_df$MemoryCD8.TcellsCD45RA._EventCounts/abundance_df$CD8.Tcell_EventCounts

naive_CD4 <- abundance_df$NaiveCD4.TCellsCD45RA._EventCounts/abundance_df$CD4.TCell_EventCounts
memory_CD4 <- abundance_df$Memory.CD4..TcellsCD45RA._EventCounts/abundance_df$CD4.TCell_EventCounts

T_subsets <- data.frame(SampleType = abundance_df$SampleType, Pair = abundance_df$Pair, 
                        naive_CD8, 
                        memory_CD8,
                        naive_CD4,
                        memory_CD4)

df <- melt(T_subsets)
ggplot(df, aes(variable, value)) + geom_boxplot(aes(colour = SampleType)) + ggtitle("T cell Abundance differences (% CD4 or CD8 T cells)") + coord_flip()
# log scale
ggplot(df, aes(variable, log10(value))) + geom_boxplot(aes(colour = SampleType)) + ggtitle("T cell Abundance differences (% CD4 or CD8 T cells), log") + coord_flip() + scale_color_brewer(palette = "Dark2") + theme_bw()
ggsave(filename = paste(figure_directory, "T cell Abundance differences.pdf", sep = ""), width = 10, height = 10,  useDingbats=FALSE)

M_df <- dplyr::filter(T_subsets, SampleType == "M")
B_df <- dplyr::filter(T_subsets, SampleType == "B")

pvals <- c()
for (i in 3:6){
  test <- t.test(M_df[,i], B_df[,i], paired = TRUE)
  pvals <- c(pvals, test$p.value)
}


# Tbet
ggplot(df, aes(variable, value)) + geom_boxplot(aes(colour = SampleType)) + ggtitle("Tbet+ CD4 T cell Abundance differences") + coord_flip() + scale_color_brewer(palette = "Dark2") + theme_bw()
ggsave(filename = paste(figure_directory, "Tbet Abundance differences.pdf", sep = ""), width = 10, height = 10,  useDingbats=FALSE)
ggplot(df, aes(variable, log10(value))) + geom_boxplot(aes(colour = SampleType)) + ggtitle("Tbet+ CD4 T cell Abundance differences, log") + coord_flip() + scale_color_brewer(palette = "Dark2") + theme_bw()
ggsave(filename = paste(figure_directory, "Tbet Abundance differences logscale.pdf", sep = ""), width = 10, height = 10,  useDingbats=FALSE)

