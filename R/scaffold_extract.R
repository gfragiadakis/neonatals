data <- scaffold:::my_load("~/Documents/neonatal/scaffold/200_clusters_new/scaffold_for_paper/030715_UC_Neonatal_UC14_M_Unstim_M_unstim_Singlets.fcs.clustered.txt.scaffold")
test <- data$graphs[["030715_UC_Neonatal_UC14_M_Unstim_M_unstim_Singlets.fcs.clustered.txt"]]

SAM_file <- read.csv("~/Documents/neonatal/scaffold/200_clusters_new/SAM_stuff/SAM_hits_thresh.csv")

V(data$graphs[["030715_UC_Neonatal_UC14_M_Unstim_M_unstim_Singlets.fcs.clustered.txt"]])[type == 2]$SAMhits <- SAM_file$hits

scaffold:::my_save(obj = data, "~/Desktop/test_file.scaffold")

data2 <- scaffold:::my_load("~/Documents/neonatal/scaffold/200_clusters_new/030715_UC_Neonatal_UC14_M_Unstim_M_unstim_Singlets.fcs.clustered.txt.scaffold")
test2 <- data$graphs[["030715_UC_Neonatal_UC14_M_Unstim_M_unstim_Singlets.fcs.clustered.txt"]]
