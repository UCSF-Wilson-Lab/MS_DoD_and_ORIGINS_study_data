#!/usr/bin/env Rscript


library(PairSeq)
library(data.table)
library(stringr)
library(dplyr)
library(janitor)
library(Biostrings)

# PATH to our local project directory is not provided
setwd("[PATH to project folder]")


# Pre-Processing ----

# 1a. Load RPK data frame
dod_rpk  <- read.csv("rpk_tables/dod_rpk.csv",stringsAsFactors = F) 
epic_rpk <- read.csv("rpk_tables/epic_rpk.csv",stringsAsFactors = F) 
gfap_rpk <- read.csv("rpk_tables/gfap_rpk.csv",stringsAsFactors = F) 
pos_rpk  <- read.csv("rpk_tables/pos_rpk.csv",stringsAsFactors = F) 
pd_rpk   <- read.csv("rpk_tables/pd_rpk.csv",stringsAsFactors = F) 


## Adding in the re sequenced samples 
re_dod_rpk          <- read.csv("rerun/dod_rpk.csv",stringsAsFactors = F) 
colnames(re_dod_rpk)<-gsub("_cgorigins","",colnames(re_dod_rpk))

missing <- setdiff(colnames(re_dod_rpk), colnames(dod_rpk))
dupes   <- intersect(colnames(re_dod_rpk), colnames(dod_rpk))
add_in  <- re_dod_rpk[,c("peptide", missing)]

## add in missing samples 
dod_rpk <- full_join(dod_rpk, add_in, by = "peptide")


# Peptide gene name table
peptide_gene_map <- read.csv("metadata/peptide_gene_mapping.csv",stringsAsFactors = F) %>% clean_names()
colnames(peptide_gene_map)[1] <- "peptide"

# merged peptide map
merged_peptide_map <- read.csv("ref_oligos/merged_peptide_map.csv", stringsAsFactors = F)
merged_peptide_map <- merged_peptide_map %>% mutate(index = row_number())

# Load meta data
metadata  <- read.csv("metadata/list_with_plate_locations.csv",stringsAsFactors = F) %>% clean_names()
prepost   <- read.csv("metadata/Case_Control_Serum.csv", stringsAsFactors = F) %>% clean_names()
metadata$type <- str_replace_all(metadata$type,"\\ ","_")
still_missing <- setdiff(metadata$serum_id, colnames(dod_rpk))


# 1b separate pre/post metadata 
collection <- c()
ser_yr <- c()

for (i in 1:nrow(metadata)){
  serID <- metadata$serum_id[i]
  if (serID %in% prepost$pre_serumid) {
    collection[i] <- "pre"
    ser_yr[i] <- prepost$pre_ser_yr[which(prepost$pre_serumid == serID)]
  } else if (serID %in% prepost$post_serumid) {
    collection[i] <- "post"
    ser_yr[i] <- prepost$post_ser_yr[which(prepost$post_serumid == serID)]    
  } else {
    collection[i] <- "NA"
    ser_yr[i]     <- "NA"
  }
}

metadata         <- cbind(metadata, collection, ser_yr)
metadata$ser_age <- metadata$ser_yr-metadata$year_birth


# 2. Format Columns/subset data
col_names <- names(dod_rpk)
col_names <- str_replace_all(col_names,"\\.","_")

cols_to_keep <- c()   
for (i in 2:length(col_names)) {
  sample    <- col_names[i]
  sample_id <- tstrsplit(sample,"_")[[1]] 
  metadata_sub <- metadata[metadata$serum_id %in% sample_id,]
  if(nrow(metadata_sub) == 0){next}
  
  type <- metadata_sub$type
  collection <- as.character(metadata_sub$collection)
  studyid <- metadata_sub$studyid
  sample_fmt <- paste(c(studyid, type, collection), collapse = "_")
  col_names[i] <- sample_fmt     
  
  cols_to_keep <- c(cols_to_keep,sample_fmt)
}

names(dod_rpk) <- col_names



# PairSeq ----
beadnames <- names(pd_rpk)
bead_cols <- beadnames[grep("bead", beadnames)]
bead_rpk  <- pd_rpk[, c("peptide", bead_cols)]

# 1. Merge both data sets (by peptide column)
df <- full_join(dod_rpk, epic_rpk, by = "peptide")
df <- full_join(df, pos_rpk, by = "peptide")
df <- full_join(df, gfap_rpk, by = "peptide")
df <- full_join(df, bead_rpk, by = "peptide")
df <- left_join(df, merged_peptide_map, by = "peptide")

# remove input dataframes
pos_cols <- names(pos_rpk)[2:6]
rm(dod_rpk, epic_rpk, gfap_rpk, pos_rpk, pd_rpk, bead_rpk)

# 2. Create unique identifier column
df$peptide_id <- paste(df$gene, df$index,sep = "_")

# Replace NA values with 0
df[is.na(df)] <- 0
any(is.na(df))


# 3. Calculate Mean RPK

# get names 
all_names         <- names(df)
gfap_cols         <- all_names[grep("gfap",all_names)]
healthy_cols      <- all_names[grep("healthy",all_names)]
healthy_pre_cols  <- healthy_cols[grep("pre", healthy_cols)]
healthy_post_cols <- healthy_cols[grep("post", healthy_cols)]  
case_pre_cols     <- all_names[grep("Case_pre",all_names)]
case_post_cols    <- all_names[grep("Case_post",all_names)]
tbi_cols          <- all_names[grep("TBI",all_names)]
mig_cols          <- all_names[grep("Migraine",all_names)]
epic_csf_cols     <- all_names[grep("_csf",all_names)]
epic_serum_cols   <- all_names[grep("_serum",all_names)]
bead_cols         <- all_names[grep("bead",all_names)]

# Mean RPK and replace 0's with median value
df$mean_ref     <- apply(df[,c(healthy_cols)], 1, mean)     # Control group
hc_avg <- median(df$mean_ref)
df$mean_ref[df$mean_ref == 0] <- hc_avg

df$mean_disease <- apply(df[,c(case_post_cols)], 1, mean)   # Cases
case_avg <- median(df$mean_disease)
df$mean_disease[df$mean_disease == 0] <- case_avg

df$mean_beads   <- apply(df[,c(bead_cols)], 1, mean)        # Beads
bead_avg <- mean(df$mean_beads)
df$mean_beads[df$mean_beads == 0] <- bead_avg

df$sd_ref <-  apply(df[,c(healthy_cols)], 1, sd) 


# 4. Calculate Fold Change
df_fc_ref    <- makeFCdf(df,
                         target_columns = healthy_cols,
                         mean_column    = "mean_ref",  onlyFC = T)
df_fc_disease <- makeFCdf(df,
                          target_columns = c(case_post_cols, case_pre_cols, tbi_cols, mig_cols, epic_csf_cols, epic_serum_cols, gfap_cols, pos_cols, bead_cols),
                          mean_column    = "mean_ref", onlyFC = T)

df_fc <- cbind(df_fc_ref,df_fc_disease)
rm(df_fc_disease, df_fc_ref)


# 5. Calculate z score
fc_names <- names(df_fc)

gfap_fc_cols <- fc_names[grep("gfap",fc_names)]
gfap_z_cols <- str_replace_all(gfap_fc_cols, "_FC","_Z")

healthy_fc_cols <- fc_names[grep("healthy",fc_names)]
healthy_z_cols <- str_replace_all(healthy_fc_cols, "_FC","_Z")

healthy_pre_fc_cols <- healthy_fc_cols[grep("pre",healthy_fc_cols)]
healthy_pre_z_cols <- str_replace_all(healthy_pre_fc_cols, "_FC","_Z")

healthy_post_fc_cols <- healthy_fc_cols[grep("post",healthy_fc_cols)]
healthy_post_z_cols <- str_replace_all(healthy_post_fc_cols, "_FC","_Z")

case_pre_fc_cols <- fc_names[grep("Case_pre",fc_names)]
case_pre_z_cols <- str_replace_all(case_pre_fc_cols, "_FC","_Z")

case_post_fc_cols <- fc_names[grep("Case_post",fc_names)]
case_post_z_cols <- str_replace_all(case_post_fc_cols, "_FC","_Z")

tbi_fc_cols <- fc_names[grep("TBI",fc_names)]
tbi_z_cols <- str_replace_all(tbi_fc_cols, "_FC","_Z")

mig_fc_cols <- fc_names[grep("Migraine",fc_names)]
mig_z_cols <- str_replace_all(mig_fc_cols, "_FC","_Z")

epic_fc_csf_cols <- fc_names[grep("_csf",fc_names)]
epic_z_csf_cols <- str_replace_all(epic_fc_csf_cols, "_FC","_Z")

epic_fc_serum_cols <- fc_names[grep("_serum",fc_names)]
epic_z_serum_cols <- str_replace_all(epic_fc_serum_cols, "_FC","_Z")

pos_fc_cols <- paste0(pos_cols, "_FC")
pos_z_cols <- str_replace_all(pos_fc_cols, "_FC","_Z")

bead_fc_cols <- paste0(bead_cols, "_FC")
bead_z_cols <- str_replace_all(pos_fc_cols, "_FC","_Z")

df_z        <- calcZscores(df_fc, target.samples = c(case_post_fc_cols, case_pre_fc_cols, tbi_fc_cols, mig_fc_cols, epic_fc_csf_cols, epic_fc_serum_cols, gfap_fc_cols, pos_fc_cols, healthy_fc_cols, bead_fc_cols), reference.samples = healthy_fc_cols)
names(df_z) <- str_replace_all(names(df_z), "_FC","_Z")
df_fc_z     <- cbind(df, df_fc, df_z)
rm(df, df_fc, df_z)







# 5. Full Parse of Candidate peptides
MIN_RPK         = 1
ZSCORE_THRESH1  = 3
ZSCORE_THRESH2  = 10
FC_THRESH1      = 10
FC_THRESH2      = 100
SUM_RPK_THRESH  = 50

start.time <- Sys.time()
cat(paste0("\n\n >> Full Parse\n\n"))
candidate_df <- fullParse(df = df_fc_z, list_of_samples = case_post_cols,
                          MIN_RPK = MIN_RPK,ZSCORE_THRESH1 = ZSCORE_THRESH1, 
                          ZSCORE_THRESH2 = ZSCORE_THRESH2,
                          FC_THRESH1 = FC_THRESH1,FC_THRESH2 = FC_THRESH2,
                          SUM_RPK_THRESH = SUM_RPK_THRESH)


end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

obj_fh <- "pairseq/fp_casepost.RData" 
save(candidate_df,file = obj_fh)

# 6. Kmer Overlap
start.time <- Sys.time()
cat(paste0("\n\n >> Kmer Overlap\n\n"))
candidate_df_fin <- findKmerOverlap.old(candidate_df,KMER_SIZE = 7)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


obj_fh <- "pairseq/kmer_casepost.RData" 
save(candidate_df_fin,file = obj_fh)


cat(paste0("\n\n>>> DONE <<<\n\n"))
