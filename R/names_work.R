old_names <- colnames(raw_df)[!colnames(raw_df) %in% annotation_columns]

names_key <- as.list(old_names)
names(names_key) <- old_names

names_key <- gsub(pattern = "Medians_", replacement = "", names_key)
names_key <- gsub(pattern = "_.*Di...", replacement = "_", names_key)
names_key <- gsub(pattern = "ERK1.2", replacement = "ERK", names_key)


# population names
NKbrights <- grep(pattern = "..1", names_key, fixed = TRUE)
names_key[NKbrights] <- gsub(pattern = "NKCD56.CD16.", replacement = "CD56hi_NK_cells", names_key[NKbrights])
names_key[NKbrights] <- gsub(pattern = "..1", replacement = "", names_key[NKbrights], fixed = TRUE)
names_key <- gsub(pattern = "NKCD56.CD16.", replacement = "CD56lo_NK_cells", names_key)

names_key <- gsub(pattern = "NaiveCD8.TcellsCD45RA.", replacement = "Naive_CD8_T_cells", names_key)
names_key <- gsub(pattern = "MemoryCD8.TcellsCD45RA.", replacement = "Memory_CD8_T_cells", names_key)
names_key <- gsub(pattern = "NaiveCD4.TCellsCD45RA.", replacement = "Naive_CD4_T_cells", names_key)
names_key <- gsub(pattern = "Memory.CD4..TcellsCD45RA.", replacement = "Memory_CD4_T_cells", names_key)
names_key <- gsub(pattern = "CD8.Tcell", replacement = "CD8_T_cells", names_key)
names_key <- gsub(pattern = "CD4.TCell", replacement = "CD4_T_cells", names_key)

names_key <- gsub(pattern = "Th1Tbet.CD4.", replacement = "Th1_T_cells", names_key)
names_key <- gsub(pattern = "gdTcells", replacement = "gd_T_cells", names_key)

names_key <- gsub(pattern = "BCellsCD19.IgM.CD3.", replacement = "B_cells", names_key)

names_key <- gsub(pattern = "ClassicalMonocytesCD14.CD16.", replacement = "Classical_monocytes", names_key)
names_key <- gsub(pattern = "Non.ClassicalMonocytesCD16.CD14.", replacement = "Nonclassical_monocytes", names_key)
names_key <- gsub(pattern = "IntermediateMonocytesCD14.CD16.", replacement = "Intermediate_monocytes", names_key)
names_key <- gsub(pattern = "M.MDSCs", replacement = "MDSCs", names_key)

names_key <- gsub(pattern = "cDCCD11c.HLA.DRhigh", replacement = "cDCs", names_key)
names_key <- gsub(pattern = "PlasmacytoidDCsHLA.DR.CD123.", replacement = "pDCs", names_key)

# remove . at the end
names_key <- gsub(pattern = ".", replacement = "", names_key, fixed = TRUE)
names(names_key) <- old_names


## replacing the names in raw_df:
replace_names <- function(names_key, df){
for (i in names(names_key)){
  if (i %in% colnames(df)){
    names(df)[names(df) == i] <- names_key[[i]]
  }
}
return(df)
}

new_df <- replace_names(names_key, raw_df)


#### 

old_names <- colnames(percents_df)[!colnames(percents_df) %in% annotation_columns]

names_key <- as.list(old_names)
names(names_key) <- old_names


# population names
NKbrights <- grep(pattern = "EventCounts.1", names_key, fixed = TRUE)
names_key[NKbrights] <- gsub(pattern = "NKCD56.CD16.", replacement = "CD56hi_NK_cells", names_key[NKbrights])
names_key[NKbrights] <- gsub(pattern = ".1", replacement = "", names_key[NKbrights], fixed = TRUE)
names_key <- gsub(pattern = "NKCD56.CD16.", replacement = "CD56lo_NK_cells", names_key)

names_key <- gsub(pattern = "NaiveCD8.TcellsCD45RA.", replacement = "Naive_CD8_T_cells", names_key)
names_key <- gsub(pattern = "MemoryCD8.TcellsCD45RA.", replacement = "Memory_CD8_T_cells", names_key)
names_key <- gsub(pattern = "NaiveCD4.TCellsCD45RA.", replacement = "Naive_CD4_T_cells", names_key)
names_key <- gsub(pattern = "Memory.CD4..TcellsCD45RA.", replacement = "Memory_CD4_T_cells", names_key)
names_key <- gsub(pattern = "CD8.Tcell", replacement = "CD8_T_cells", names_key)
names_key <- gsub(pattern = "CD4.TCell", replacement = "CD4_T_cells", names_key)

names_key <- gsub(pattern = "Th1Tbet.CD4.", replacement = "Th1_T_cells", names_key)
names_key <- gsub(pattern = "gdTcells", replacement = "gd_T_cells", names_key)

names_key <- gsub(pattern = "BCellsCD19.IgM.CD3.", replacement = "B_cells", names_key)

names_key <- gsub(pattern = "ClassicalMonocytesCD14.CD16.", replacement = "Classical_monocytes", names_key)
names_key <- gsub(pattern = "Non.ClassicalMonocytesCD16.CD14.", replacement = "Nonclassical_monocytes", names_key)
names_key <- gsub(pattern = "IntermediateMonocytesCD14.CD16.", replacement = "Intermediate_monocytes", names_key)
names_key <- gsub(pattern = "M.MDSCs", replacement = "MDSCs", names_key)

names_key <- gsub(pattern = "cDCCD11c.HLA.DRhigh", replacement = "cDCs", names_key)
names_key <- gsub(pattern = "PlasmacytoidDCsHLA.DR.CD123.", replacement = "pDCs", names_key)

names_key <- gsub(pattern = "EventCounts", replacement = "Percent_Total", names_key)

names(names_key) <- old_names
new_percent_df <- replace_names(names_key, percents_df)
new_percent <- replace_names(names_key, percents)

