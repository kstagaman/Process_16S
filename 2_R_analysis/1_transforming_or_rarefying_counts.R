# 1_transforming_or_rarefying_counts.R

workdir <- "~/Dropbox/Rag1_manuscript/Stagaman_etal_ISME_materials/Supplemental_code/2_R_analysis/"
setwd(workdir)

library(DESeq2)
library(phangorn)
library(phyloseq)

dirs <- c("B_Plots", 
          "C_Text_tables",
          "D_LefSe_files",
          "E_Saved_objects")
for (dir in dirs) {
  if (dir.exists(dir) == FALSE) {
    dir.create(dir)
  }
}
rare_seed=874869
map_file <- "A_Input_files/current_mapping_file_24Sep2015.txt"
otu_file <- "A_Input_files/all_smpls.otu_tbl.txt"
tree_file <- "A_Input_files/all_smpls.otus.tree"
tax_file <- "A_Input_files/all_smpls.taxa_tbl.phyloseq.txt"

save_nofood_rawCts_file <- "E_Saved_objects/rawCts_nofood"
save_nofood_rarCts_file <- "E_Saved_objects/rarCts_nofood"
save_gut_rarCts_file <- "E_Saved_objects/rarCts_gut"
save_water_rarCts_file <- "E_Saved_objects/rarCts_water"
save_nofood_varStab_file <- "E_Saved_objects/varStab_nofood"

############

# Import files and create phyloseq object
otu_raw <- read.table(file=otu_file, sep = "\t", row.names = 1, header = TRUE, comment.char="")
otus <- otu_table(otu_raw, taxa_are_rows=FALSE)

map_raw <- read.table(file=map_file, sep = "\t", row.names = 1, header = TRUE, comment.char="")
map_raw$dpf <- gsub("9", "09", map_raw$dpf)
map_raw$dpf <- as.factor(map_raw$dpf)
map <- sample_data(map_raw)

tree <- read_tree(tree_file)
tree.mp <- midpoint(tree)

tax_raw <- read.table(file=tax_file, sep="\t",row.names=1, header=TRUE)
tax_raw <- as.matrix(tax_raw)
taxa <- tax_table(tax_raw)

rag1.physeq <- phyloseq(otus, map, tree.mp, taxa)

# Remove taxa matching to chloroplast (probably food)
no.chloro.physeq <- subset_taxa(rag1.physeq, Class != "Chloroplast" & Phylum != "Chloroplast")

############

no.food.physeq <- subset_samples(no.chloro.physeq, smpl.type != "food")
no.food.physeq <- subset_samples(no.food.physeq, host.gt != "ht")
no.food.physeq <- subset_samples(no.food.physeq, sample_sums(no.food.physeq) > 0)

gut.physeq <- subset_samples(no.food.physeq, smpl.type=="gut")
gut.physeq <- subset_samples(gut.physeq, sample_sums(gut.physeq) > 10000)
water.physeq <- subset_samples(no.food.physeq, smpl.type=="water")

no.food.physeq <- merge_phyloseq(gut.physeq, water.physeq)
no.food.physeq <- prune_taxa(taxa_sums(no.food.physeq) > 0, no.food.physeq)
sample_data(no.food.physeq)$housing <- factor(sample_data(no.food.physeq)$housing, 
                                       levels=c("isolated", "cohoused"))
sample_data(no.food.physeq)$host.gt <- factor(sample_data(no.food.physeq)$host.gt, 
                                       levels=c("wt", "ko", "none"))
sample_data(no.food.physeq)$tank.gt <- factor(sample_data(no.food.physeq)$tank.gt, 
                                       levels=c("all.wt", "mix", "all.ko"))
save(no.food.physeq, file=save_nofood_rawCts_file)

gut.physeq <- rarefy_even_depth(gut.physeq, replace=FALSE, rngseed=rare_seed)
save(gut.physeq, file=save_gut_rarCts_file)
water.physeq <- rarefy_even_depth(water.physeq, replace=FALSE, rngseed=rare_seed, trimOTUs=FALSE)
save(water.physeq, file=save_water_rarCts_file)

############
# Create phyloseq object with variance-stabilized OTU counts

rm(no.food.physeq)
load(save_nofood_rawCts_file, verbose=TRUE)

no.food.physeq.dds <- phyloseq_to_deseq2(no.food.physeq, ~ dpf + housing + host.gt)
no.food.physeq.dds <- estimateSizeFactors(no.food.physeq.dds)
no.food.physeq.dds <- estimateDispersions(no.food.physeq.dds)
no.food.physeq.vst <- getVarianceStabilizedData(no.food.physeq.dds)
otu_table(no.food.physeq) <- otu_table(no.food.physeq.vst, taxa_are_rows = TRUE)
otu_table(no.food.physeq)[otu_table(no.food.physeq) < 0] <- 0
save(no.food.physeq, file=save_nofood_varStab_file)