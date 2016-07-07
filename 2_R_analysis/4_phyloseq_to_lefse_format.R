# 4_phyloseq_to_lefse_format.R

# workdir <- ""
# setwd(workdir)

library(phyloseq)
source("0_Support_scripts/phyloseq_to_lefse.R")

save_nofood_varStab_file <- "E_Saved_objects/varStab_nofood"

############

load(save_nofood_varStab_file, verbose=TRUE)

gut.d75.iso.physeq <- subset_samples(no.food.physeq, smpl.type=="gut" & dpf=="75" & housing=="isolated")
gut.d75.iso.physeq <- prune_taxa(taxa_sums(gut.d75.iso.physeq) > 0, gut.d75.iso.physeq)

gut.d75.coh.physeq <- subset_samples(no.food.physeq, smpl.type=="gut" & dpf=="75" & housing=="cohoused")
gut.d75.coh.physeq <- prune_taxa(taxa_sums(gut.d75.coh.physeq) > 0, gut.d75.coh.physeq)

phyloseq_to_lefse(phyloseq.object=gut.d75.iso.physeq, 
                  row_names=c("host.gt", "smpl"), 
                  name="D_LefSe_files/gut_d75_seg")
phyloseq_to_lefse(phyloseq.object=gut.d75.coh.physeq, 
                  row_names=c("host.gt", "smpl"), 
                  name="D_LefSe_files/gut_d75_mix")