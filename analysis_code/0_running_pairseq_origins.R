#!/usr/bin/env Rscript

library(tidyverse)
library(ggpubr)
library(PairSeq)
library(data.table)
library(stringr)
library(dplyr)
library(janitor)
library(Biostrings)
library(pheatmap)
library(wesanderson)
library(treemapify)
library(ComplexHeatmap)
library(RColorBrewer)
library(viridis)

# PATH to our local project directory is not provided
setwd("[PATH to project folder]")


# 1a. Load MEEBO RPK data frame
bead_rpk    <- read.csv("data_origins/beads_rpk.csv",stringsAsFactors = F) 
nid_rpk     <- read.csv("data_origins/nid_rpk.csv",stringsAsFactors = F) 
nybb_rpk    <- read.csv("data_origins/nybb_rpk.csv",stringsAsFactors = F) 
origins_rpk <- read.csv("data_origins/origins_rpk.csv",stringsAsFactors = F) 
gfap_rpk    <- read.csv("data_origins/gfap_rpk.csv",stringsAsFactors = F) 

## 2. Load and metadata
metadata <- read.csv("data_origins/ORIGINSEPICPatientSamples.csv",stringsAsFactors = F) %>% clean_names()
colnames(metadata)[1] <- "sample"

ond_ids  <- read.csv("data_origins/OND_ID_list.csv",stringsAsFactors = F)
ond      <- ond_ids$ID
ond      <- unique(ond)

# Get Sample Vectors
count_df <- metadata %>% filter(final_diagnosis == "MS_csf") 
ms_csf_cols <- count_df[, "sample"]

count_df <- metadata %>% filter(final_diagnosis == "MS_serum") 
ms_serum_cols <- count_df[, "sample"]

count_df <- metadata %>% filter(final_diagnosis == "OND_csf") 
ond_csf_cols <- count_df[, "sample"]

count_df <- metadata %>% filter(final_diagnosis == "OND_serum") 
ond_serum_cols <- count_df[, "sample"]

count_df <- metadata %>% filter(final_diagnosis == "Healthy_csf") 
healthy_csf_cols <- count_df[, "sample"]

count_df <- metadata %>% filter(final_diagnosis == "Healthy_serum") 
healthy_serum_cols <- count_df[, "sample"]

count_df <- metadata %>% filter(final_diagnosis == "CIS_serum") 
CIS_serum_cols <- count_df[, "sample"]

count_df <- metadata %>% filter(final_diagnosis == "CIS_csf") 
CIS_csf_cols <- count_df[, "sample"]

count_df <- metadata %>% filter(final_diagnosis ==  "RIS_serum") 
RIS_serum_cols <- count_df[, "sample"]

count_df <- metadata %>% filter(final_diagnosis == "RIS_csf") 
RIS_csf_cols <- count_df[, "sample"]

count_df <- metadata %>% filter(final_diagnosis ==  "Unknown_serum") 
unknown_serum_cols <- count_df[, "sample"]

count_df <- metadata %>% filter(final_diagnosis == "Unknown_csf") 
unknown_csf_cols <- count_df[, "sample"]

rm(count_df)


## 3.Get all relevant sample lists
nid_cols  <- colnames(nid_rpk)[-1]
bead_cols <- colnames(bead_rpk)[-1]
nybb_cols <- colnames(nybb_rpk)[-1]
origins_cols <- colnames(origins_rpk)[-1]

#MS
origins_ms_csf_cols   <- paste0("X", ms_csf_cols, "_csf_cgorigins")
origins_ms_serum_cols <- paste0("X", ms_serum_cols, "_serum_cgorigins")
origins_ms_csf_cols   <- intersect(origins_cols, origins_ms_csf_cols)
origins_ms_serum_cols <- intersect(origins_cols, origins_ms_serum_cols)
origins_ms_cols       <- c(origins_ms_csf_cols, origins_ms_serum_cols)

#healthy controls
origins_healthy_csf_cols   <- c(paste0("X", healthy_csf_cols, "_csf_cgorigins"))
origins_healthy_serum_cols <- paste0("X", healthy_serum_cols, "_serum_cgorigins")
origins_healthy_csf_cols   <- intersect(origins_cols, origins_healthy_csf_cols)
origins_healthy_serum_cols <- intersect(origins_cols, origins_healthy_serum_cols)


#OND
origins_ond_csf_cols   <- c(paste0("X", ond_csf_cols, "_csf_cgorigins"))
origins_ond_serum_cols <- paste0("X", ond_serum_cols, "_serum_cgorigins")
origins_ond_csf_cols   <- intersect(origins_cols, origins_ond_csf_cols)
origins_ond_serum_cols <- intersect(origins_cols, origins_ond_serum_cols)

#CIS
origins_CIS_csf_cols   <- paste0("X", CIS_csf_cols, "_csf_cgorigins")
origins_CIS_serum_cols <- paste0("X", CIS_serum_cols, "_serum_cgorigins")
origins_CIS_csf_cols   <- intersect(origins_cols, origins_CIS_csf_cols)
origins_CIS_serum_cols <- intersect(origins_cols, origins_CIS_serum_cols)

origins_RIS_csf_cols   <- paste0("X", RIS_csf_cols, "_csf_cgorigins")
origins_RIS_serum_cols <- paste0("X", RIS_serum_cols, "_serum_cgorigins")
origins_RIS_csf_cols   <- intersect(origins_cols, origins_RIS_csf_cols)
origins_RIS_serum_cols <- intersect(origins_cols, origins_RIS_serum_cols)

origins_unknown_csf_cols   <- paste0("X", unknown_csf_cols, "_csf_cgorigins")
origins_unknown_serum_cols <- paste0("X", unknown_serum_cols, "_serum_cgorigins")
origins_unknown_csf_cols   <- intersect(origins_cols, origins_unknown_csf_cols)
origins_unknown_serum_cols <- intersect(origins_cols, origins_unknown_serum_cols)

all_origins <- c(origins_healthy_csf_cols, origins_healthy_serum_cols, origins_ond_csf_cols, origins_ond_serum_cols, origins_ms_csf_cols, origins_ms_serum_cols,
                 origins_CIS_csf_cols , origins_CIS_serum_cols, origins_RIS_csf_cols, origins_RIS_serum_cols,
                 origins_unknown_csf_cols,  origins_unknown_serum_cols)
mycraw <- setdiff(origins_cols, all_origins)

nid_csf_cols    <- nid_cols[grep("_csf_",nid_cols)]
ctrl_csf_cols   <- c(origins_healthy_csf_cols, origins_ond_csf_cols, nid_csf_cols )
ctrl_serum_cols <- c(nybb_cols, origins_healthy_serum_cols, origins_ond_serum_cols)
gfap_cols       <- c("GFAP_1_cgorigins", "GFAP_2_cgorigins", "GFAP_3_cgorigins", "GFAP_4_cgorigins")


# Peptide gene name table
peptide_gene_map <- read.csv("metadata/peptide_gene_mapping.csv",stringsAsFactors = F) %>% clean_names()
colnames(peptide_gene_map)[1] <- "peptide"

#merged peptide map
merged_peptide_map <- read.csv("ref_oligos/merged_peptide_map.csv", stringsAsFactors = F)
merged_peptide_map <- merged_peptide_map %>% mutate(index = row_number())

# Load meta data
nid_rpk <- nid_rpk[, c("peptide", nid_csf_cols)]

# 1. Merge both data sets (by peptide column)
df <- full_join(origins_rpk, nid_rpk, by = "peptide")
df <- full_join(df, nybb_rpk, by = "peptide")
df <- full_join(df, bead_rpk, by = "peptide")
df <- full_join(df, gfap_rpk, by = "peptide")
df <- left_join(df, merged_peptide_map, by = "peptide")

#remove input dataframes
rm(origins_rpk, nid_rpk, bead_rpk, nybb_rpk, gfap_rpk)

# 2. Create unique identifier column
df$peptide_id <- paste(df$gene, df$index,sep = "_")

#Replace NA values with 0
df[is.na(df)] <- 0
any(is.na(df))


# 3. Calculate Mean RPK

# Mean RPK and replace 0's with median value

##CSF
df$mean_csf_ref     <- apply(df[,c(ctrl_csf_cols)], 1, mean)     # Control group
hc_avg <- median(df$mean_csf_ref)
df$mean_csf_ref[df$mean_csf_ref == 0] <- hc_avg

df$mean_csf_disease <- apply(df[,c(origins_ms_csf_cols)], 1, mean)   # Cases
case_avg <- median(df$mean_csf_disease)
df$mean_csf_disease[df$mean_csf_disease == 0] <- case_avg

##Serum
df$mean_serum_ref     <- apply(df[,c(ctrl_serum_cols)], 1, mean)     # Control group
hc_avg <- median(df$mean_serum_ref)
df$mean_serum_ref[df$mean_serum_ref == 0] <- hc_avg

df$mean_serum_disease <- apply(df[,c(origins_ms_serum_cols)], 1, mean)   # Cases
case_avg <- median(df$mean_serum_disease)
df$mean_serum_disease[df$mean_serum_disease == 0] <- case_avg

##Beads
df$mean_beads   <- apply(df[,c(bead_cols)], 1, mean)        # Beads
bead_avg <- mean(df$mean_beads)
df$mean_beads[df$mean_beads == 0] <- bead_avg


# 4. Calculate Fold Change

#CSF
df_fc_csf_ref    <- makeFCdf(df,
                             target_columns = ctrl_csf_cols,
                             mean_column    = "mean_csf_ref",  onlyFC = T)
df_fc_csf_disease <- makeFCdf(df,
                              target_columns = c(origins_ms_csf_cols),
                              mean_column    = "mean_csf_ref", onlyFC = T)

#Serum
df_fc_serum_ref    <- makeFCdf(df,
                               target_columns = ctrl_serum_cols,
                               mean_column    = "mean_serum_ref",  onlyFC = T)
df_fc_serum_disease <- makeFCdf(df,
                                target_columns = c(origins_ms_serum_cols, bead_cols),
                                mean_column    = "mean_serum_ref", onlyFC = T)


#Bind and delete
df_fc <- cbind(df_fc_csf_ref, df_fc_csf_disease, df_fc_serum_ref, df_fc_serum_disease)


rm(df_fc_csf_ref,df_fc_csf_disease, df_fc_serum_ref, df_fc_serum_disease)


# 5. Calculate z score
ctrl_fc_csf_cols <- paste0(ctrl_csf_cols, "_FC")
ctrl_z_csf_cols <- str_replace_all(ctrl_fc_csf_cols, "_FC","_Z")

ctrl_fc_serum_cols <- paste0(ctrl_serum_cols, "_FC")
ctrl_z_serum_cols <- str_replace_all(ctrl_fc_serum_cols, "_FC","_Z")

origins_ms_fc_csf_cols <- paste0(origins_ms_csf_cols, "_FC")
origins_ms_z_csf_cols <- str_replace_all(origins_ms_fc_csf_cols, "_FC","_Z")

origins_ms_fc_serum_cols <- paste0(origins_ms_serum_cols, "_FC")
origins_ms_z_serum_cols <- str_replace_all(origins_ms_fc_serum_cols, "_FC","_Z")

bead_fc_cols <- paste0(bead_cols, "_FC")
bead_z_cols <- str_replace_all(bead_fc_cols, "_FC","_Z")

#csf
df_csf_z <- calcZscores(df_fc, target.samples = c(origins_ms_fc_csf_cols, ctrl_fc_csf_cols), reference.samples = ctrl_fc_csf_cols)
names(df_csf_z) <- str_replace_all(names(df_csf_z), "_FC","_Z")

#serum
df_serum_z <- calcZscores(df_fc, target.samples = c(origins_ms_fc_serum_cols, ctrl_fc_serum_cols, bead_fc_cols), reference.samples = ctrl_fc_serum_cols)
names(df_serum_z) <- str_replace_all(names(df_serum_z), "_FC","_Z")

#Bind and delete
df_z <- cbind(df_csf_z, df_serum_z)
rm(df_csf_z, df_serum_z)

#Stick the petide info back on
df_fc_z <- cbind(df, df_fc, df_z)

#remove input dataframes

rm(df, df_fc, df_z)


# 5. Full Parse of Candidate peptides

#CSF
MIN_RPK         = 1
ZSCORE_THRESH1  = 3
ZSCORE_THRESH2  = 10
FC_THRESH1      = 10
FC_THRESH2      = 100
SUM_RPK_THRESH  = 50



start.time <- Sys.time()
cat(paste0("\n\n >> Full Parse CSF\n\n"))
candidate_df_csf <- fullParse(df = df_fc_z, list_of_samples = origins_ms_csf_cols,
                          MIN_RPK = MIN_RPK, ZSCORE_THRESH1 = ZSCORE_THRESH1,
                          ZSCORE_THRESH2 = ZSCORE_THRESH2,
                          FC_THRESH1 = FC_THRESH1,FC_THRESH2 = FC_THRESH2,
                          SUM_RPK_THRESH = SUM_RPK_THRESH)


end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

obj_fh <- "pairseq/origins_csf_fp.RData" 
save(candidate_df_csf,file = obj_fh)

# 6. CSF Kmer Overlap
start.time <- Sys.time()
cat(paste0("\n\n >> Kmer Overlap CSF\n\n"))
candidate_df_csf_fin <- findKmerOverlap.old(candidate_df_csf,KMER_SIZE = 7)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


obj_fh <- "pairseq/origins_csf_kmer.RData" 
save(candidate_df_csf_fin,file = obj_fh)

##     Serum
MIN_RPK         = 1
ZSCORE_THRESH1  = 3
ZSCORE_THRESH2  = 10
FC_THRESH1      = 10
FC_THRESH2      = 100
SUM_RPK_THRESH  = 50



start.time <- Sys.time()
cat(paste0("\n\n >> Full Parse Serum\n\n"))
candidate_df_serum <- fullParse(df = df_fc_z, list_of_samples = origins_ms_serum_cols,
                          MIN_RPK = MIN_RPK, ZSCORE_THRESH1 = ZSCORE_THRESH1,
                          ZSCORE_THRESH2 = ZSCORE_THRESH2,
                          FC_THRESH1 = FC_THRESH1,FC_THRESH2 = FC_THRESH2,
                          SUM_RPK_THRESH = SUM_RPK_THRESH)


end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

obj_fh <- "pairseq/origins_serum_fp.RData" 
save(candidate_df_serum,file = obj_fh)

# 6. Serum Kmer Overlap
start.time <- Sys.time()
cat(paste0("\n\n >> Kmer Overlap serum\n\n"))
candidate_df_serum_fin <- findKmerOverlap.old(candidate_df_serum,KMER_SIZE = 7)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

obj_fh <- "pairseq/origins_serum_kmer.RData" 
save(candidate_df_serum_fin,file = obj_fh)


cat(paste0("\n\n>>> DONE <<<\n\n"))
