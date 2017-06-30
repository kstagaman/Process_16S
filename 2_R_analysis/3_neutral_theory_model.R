# 3_neutral_theory_model.R

# workdir <- ""
# setwd(workdir)

library(car)
library(DTK)
library(multcompView)
library(plyr)
library(doBy)
library(gridExtra)
library(RColorBrewer)
library(qdap)
library(ggplot2)
theme_set(theme_bw(base_size=18))
library(phyloseq)
source("0_Support_scripts/generate_label_df_errorbars.R")
source("0_Support_scripts/Supplemental_Code_1.R")
source("0_Support_scripts/neutral_theory_data.R")

save_gut_rarCts_file <- "E_Saved_objects/rarCts_gut" 
partitioned_physeq_file <- "E_Saved_objects/partitioned_physeq_obj_byTank-T" 

brewer.Set1 <- brewer.pal(9, "Set1")
brewer.Paired <- brewer.pal(12, "Paired")
housing.labs <- c("Segregated", "Mixed")
housing.shps <- c(19, 4)
host.gt.labs <- c("Wild type", expression(italic(paste("rag1"))^-{}))
gt.cols <- c(brewer.Paired[2], brewer.Set1[1]) # blue, red

############

# Generate data frame for neutral model fit plot
do.rsq <- FALSE
load(save_gut_rarCts_file, verbose=TRUE)
gut.data <- as(sample_data(gut.physeq), "data.frame")
i <- 0
n.vec <- NULL
for(house.type in levels(gut.data$housing)) {
    for(genotype in levels(gut.data$host.gt)) {
        i <- i+1
        print(paste(house.type, "-", genotype, ": ", 
                    nrow(subset(gut.data, housing==house.type & host.gt==genotype)), sep=""))
        n.vec[i] <- nrow(subset(gut.data, housing==house.type & host.gt==genotype))
    }
}
min.smpls <- min(n.vec) - 2
paste("min samples:", min.smpls)
reps <- 99

allGut.rsq.file <- paste("E_Saved_objects/Rsq_fit_data_allGutSource_minSmpls", 
                         min.smpls, 
                         "_reps", 
                         reps, 
                         sep="")
pool_otu_file <- paste("E_Saved_objects/rarCts_gut_otu_tbl.txt", sep="")
write.table(otu_table(gut.physeq), 
            file=pool_otu_file, 
            sep="\t")
pool.otu.tbl <- as.matrix(read.table(pool_otu_file))

if(do.rsq == TRUE) {
    r.sq.data <- NULL
    for(house.type in unique(as.character(gut.data$housing))) {
        for(genotype in unique(as.character(gut.data$host.gt))) {
            sub.data <- subset(gut.data, housing==house.type & host.gt==genotype)
            smpl.names <- row.names(sub.data)
            for(i in c(1:reps)) {
                keep.smpls <- sample(smpl.names, min.smpls, replace=FALSE)
                curr.physeq <- prune_samples(keep.smpls, gut.physeq)
                curr.physeq <- prune_taxa(taxa_sums(curr.physeq) > 0, curr.physeq)
                otu_file <- paste("E_Saved_objects/otu_tbl_", 
                                  house.type, 
                                  "_", 
                                  genotype, 
                                  i, 
                                  ".txt", 
                                  sep="")
                tax_file <- paste("E_Saved_objects/tax_tbl_", 
                                  house.type, 
                                  "_", 
                                  genotype, 
                                  i, 
                                  ".txt", 
                                  sep="")
                write.table(otu_table(curr.physeq), file=otu_file, sep="\t")
                write.table(tax_table(curr.physeq), file=tax_file, sep="\t")
                otu.tbl <- as.matrix(read.table(otu_file))
                tax.tbl <- as.matrix(read.table(tax_file))
                
                fit <- sncm.fit(otu.tbl, pool=pool.otu.tbl, taxon=tax.tbl, stats=TRUE)
                fit.df <- cbind(smpl.type=unique(sample_data(curr.physeq)$smpl.type),
                                housing=unique(sample_data(curr.physeq)$housing),
                                host.gt=unique(sample_data(curr.physeq)$host.gt),
                                fit)
                r.sq.data <- rbind(r.sq.data, fit.df)
            }
        }
    }
    
    save(r.sq.data, file=allGut.rsq.file)
}

############

# Plot R-squared values contained in data frame generated above.
load(allGut.rsq.file, verbose=TRUE)
allGut.rsq.plot <- ggplot(r.sq.data, aes(x=housing, y=Rsqr))
allGut.rsq.geoms <- allGut.rsq.plot + 
    geom_point(aes(shape=housing, color=host.gt, fill=host.gt, alpha=housing), 
               position=position_jitterdodge(jitter.width=0.7, dodge.width=1), 
               size=3) + 
    scale_alpha_discrete(range=c(0.4, 0.7), guide=FALSE) + 
    stat_summary(aes(group=interaction(housing, host.gt)), 
                 fun.data="mean_se", 
                 geom="errorbar", 
                 position="dodge", 
                 width=1, 
                 size=1, 
                 color="gray21") + 
    scale_x_discrete(name="Housing", 
                     labels=housing.labs) + 
    scale_shape_manual(values=housing.shps) + 
    scale_linetype_manual(values=housing.ltys) + 
    scale_color_manual(name="Host genotype", 
                       values=gt.cols, 
                       labels=host.gt.labs) + 
    labs(y="Goodness of fit\n(Generalized R2)") + 
    theme(strip.background=element_blank(),
          strip.text=element_text(hjust=0),
          legend.position="right", 
          plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"),
          plot.background=element_rect(fill=NA),
          legend.text.align=0) + 
    guides(fill=FALSE, shape=FALSE, linetype=FALSE)
factors <- unique(as.character(allGut.rsq.geoms$layers[[1]]$mapping))
measure <- as.character(allGut.rsq.geoms$mapping$y)
names(factors) <- NULL
factors.eq <- paste(factors, collapse="*")
aov.eq <- paste(measure, "~", factors.eq)
sig.labels.list <- generate_label_df_errorbars(orig.data=r.sq.data, 
                                               ggplot.obj=allGut.rsq.geoms, 
                                               aov.eq)
sig.labels <- sig.labels.list$labels.df
save.allGut.rsq.plot <- allGut.rsq.geoms + 
    geom_text(data=sig.labels, size=5, aes(y=V1, x=x.coords, label=labels))
# save.allGut.rsq.plot
ggsave(save.allGut.rsq.plot,
       file=paste(gsub("E_Saved_objects", "B_Plots", allGut.rsq.file), "_plot.pdf", sep=""),
       units="mm",
       width=86,
       height=54)

############

# Generate lists of OTUs in each partition (by Tank)
source("0_Support_scripts/neutral_theory_data.R")
load(save_gut_rarCts_file, verbose=TRUE)
gut.d75.physeq <- subset_samples(gut.physeq, dpf=="75")
gut.d75.physeq <- prune_taxa(taxa_sums(gut.d75.physeq) > 0, gut.d75.physeq)

if(grepl("byTank-T", partitioned_physeq_file)) {
    neutral.theory.data <- neutral_theory_data(physeq.object=gut.d75.physeq, 
                                               by.tank=TRUE, 
                                               source.pool=TRUE)
    prop.part.eval <- TRUE
} else if (grepl("byTank-F", partitioned_physeq_file)) {
    neutral.theory.data <- neutral_theory_data(physeq.object=gut.d75.physeq, 
                                               by.tank=FALSE, 
                                               source.pool=TRUE)
    prop.part.eval <- FALSE
}

all.neutral.theory <- neutral.theory.data$all
save(all.neutral.theory, file="E_Saved_objects/all_neutral_theory")
nonNeutral.otus <- as.character(unique(neutral.theory.data$nonNeutral$otu))
neutral.otus <- as.character(unique(neutral.theory.data$neutral$otu))
neutral.otus.data <- neutral.theory.data$neutral
nonNeutral.otus.data <- neutral.theory.data$nonNeutral
treatments <- levels(neutral.otus.data$name)

neutral.otus.list <- list()
for(treatment in treatments) {
    sub.data <- subset(neutral.otus.data, name == treatment)
    sub.list <- list(as.character(sub.data$otu))
    names(sub.list) <- treatment
    write.table(sub.list, 
                file=paste("E_Saved_objects/", 
                           gsub(" ", "_", 
                                paste(treatment, 
                                      "neutral_otu_list.txt", 
                                      sep="_")), 
                           sep=""), 
                quote=FALSE, 
                sep="\t",
                row.names=FALSE)
    neutral.otus.list <- c(neutral.otus.list, sub.list)
}
underRep.otus.list <- list()
for(treatment in treatments) {
    sub.data <- subset(nonNeutral.otus.data, name == treatment & selection == "Under represented")
    sub.list <- list(as.character(sub.data$otu))
    names(sub.list) <- treatment
    write.table(sub.list, 
                file=paste("E_Saved_objects/", 
                           gsub(" ", "_", 
                                paste(treatment, 
                                      "selected_against_otu_list.txt", 
                                      sep="_")), 
                           sep=""), 
                quote=FALSE, 
                sep="\t",
                row.names=FALSE)
    underRep.otus.list <- c(underRep.otus.list, sub.list)
}
overRep.otus.list <- list()
for(treatment in treatments) {
    sub.data <- subset(nonNeutral.otus.data, name == treatment & selection == "Over represented")
    sub.list <- list(as.character(sub.data$otu))
    names(sub.list) <- treatment
    write.table(sub.list, 
                file=paste("E_Saved_objects/", 
                           gsub(" ", "_", 
                                paste(treatment, 
                                      "selected_for_otu_list.txt", 
                                      sep="_")), 
                           sep=""), 
                quote=FALSE, 
                sep="\t",
                row.names=FALSE)
    overRep.otus.list <- c(overRep.otus.list, sub.list)
}


############

# Using the lists above, I partition the community of each individual sample and create a new phyloseq object with that OTU table
curr1 <- c("Knockout", "Wild type", "Mixed", "Segregated")
simple1 <- c("ko", "wt", "cohoused", "isolated")

partitioned.physeq <- NULL
for(i in c(1:length(names(neutral.otus.list)))){
    name <- names(neutral.otus.list)[i]
    name <- mgsub(curr1, simple1, name)
    name.split <- strsplit(name, "-", fixed=TRUE)
    if (grepl("byTank-T", partitioned_physeq_file)) {
        subset.physeq <- subset_samples(gut.d75.physeq, 
                                        host.gt==name.split[[1]][1] & 
                                            housing==name.split[[1]][2] & 
                                            tank==name.split[[1]][3])
    } else if (grepl("byTank-F", partitioned_physeq_file)) {
        subset.physeq <- subset_samples(gut.d75.physeq, 
                                        host.gt==name.split[[1]][1] & 
                                            housing==name.split[[1]][2])
    }
    sub.smpl.df <- as.data.frame(sample_data(subset.physeq))
    sub.otu.df <- as.data.frame(otu_table(subset.physeq))
    
    # Over represented partition
    overRep.smpl.df <- sub.smpl.df
    overRep.smpl.df <- cbind(overRep.smpl.df, 
                             "partition"=rep("overRep", 
                                             length(overRep.smpl.df$dpf)))
    overRep.otu.tbl <- sub.otu.df
    row.names(overRep.smpl.df) <- paste(row.names(overRep.smpl.df), "-overRep", sep="")
    row.names(overRep.otu.tbl) <- paste(row.names(overRep.otu.tbl), "-overRep", sep="")
    zeroed.otus <- names(sub.otu.df)[!(names(sub.otu.df) %in% overRep.otus.list[[i]])]
    overRep.otu.tbl[, zeroed.otus] <- 0
    for.physeq <- phyloseq(sample_data(overRep.smpl.df), 
                           otu_table(overRep.otu.tbl, taxa_are_rows=FALSE), 
                           tax_table(subset.physeq), 
                           phy_tree(subset.physeq))
    partitioned.physeq <- merge_phyloseq(partitioned.physeq, for.physeq)
    
    # Under represented partition
    underRep.smpl.df <- sub.smpl.df
    underRep.smpl.df <- cbind(underRep.smpl.df, 
                              "partition"=rep("underRep", 
                                              length(underRep.smpl.df$dpf)))
    underRep.otu.tbl <- sub.otu.df
    row.names(underRep.smpl.df) <- paste(row.names(underRep.smpl.df), "-underRep", sep="")
    row.names(underRep.otu.tbl) <- paste(row.names(underRep.otu.tbl), "-underRep", sep="")
    zeroed.otus <- names(sub.otu.df)[!(names(sub.otu.df) %in% underRep.otus.list[[i]])]
    underRep.otu.tbl[, zeroed.otus] <- 0
    agn.physeq <- phyloseq(sample_data(underRep.smpl.df), 
                           otu_table(underRep.otu.tbl, taxa_are_rows=FALSE), 
                           tax_table(subset.physeq), 
                           phy_tree(subset.physeq))
    partitioned.physeq <- merge_phyloseq(partitioned.physeq, agn.physeq)
    
    # NEUTRAL partition
    neutral.smpl.df <- sub.smpl.df
    neutral.smpl.df <- cbind(neutral.smpl.df, 
                             "partition"=rep("neutral", 
                                             length(neutral.smpl.df$dpf)))
    neutral.otu.tbl <- sub.otu.df
    row.names(neutral.smpl.df) <- paste(row.names(neutral.smpl.df), "-neutral", sep="")
    row.names(neutral.otu.tbl) <- paste(row.names(neutral.otu.tbl), "-neutral", sep="")
    zeroed.otus <- names(sub.otu.df)[!(names(sub.otu.df) %in% neutral.otus.list[[i]])]
    neutral.otu.tbl[, zeroed.otus] <- 0
    neu.physeq <- phyloseq(sample_data(neutral.smpl.df), 
                           otu_table(neutral.otu.tbl, taxa_are_rows=FALSE), 
                           tax_table(subset.physeq), 
                           phy_tree(subset.physeq))
    partitioned.physeq <- merge_phyloseq(partitioned.physeq, neu.physeq)
}
partitioned.physeq <- prune_samples(sample_sums(partitioned.physeq) > 0, partitioned.physeq)
save(partitioned.physeq, file=partitioned_physeq_file)

############

### Making ordinations and running PERMANOVAS
load(partitioned_physeq_file, verbose=TRUE)
partitions <- c("neutral", "overRep", "underRep")

partOrd.list <- list()
betadisp.list <- list()
for(i in c(1:length(partitions))) {
    physeq <- subset_samples(partitioned.physeq, partition == partitions[i])
    physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq)
    ord <- phyloseq::ordinate(physeq, 
                              method="NMDS", 
                              distance=dist.method, 
                              binary=dist.binary, 
                              try.max=99)
    ord.data <- plot_ordination(physeq, ord, justDF=TRUE)
    ord.data$partition <- factor(mgsub(c("neutral", 
                                         "overRep", 
                                         "underRep"), 
                                       c("Neutral", 
                                         "Over represented", 
                                         "Under represented"),
                                       ord.data$partition),
                                 levels=c("Neutral", 
                                          "Over represented", 
                                          "Under represented"))
    ord.plot <- ggplot(ord.data, aes(x=NMDS1, y=NMDS2))
    save.ord.plot <- ord.plot +
        geom_point(aes(shape=housing, color=host.gt), size=1) + 
        scale_shape_manual(values=as.numeric(mgsub(c("isolated", "cohoused"), 
                                                   housing.shps, 
                                                   levels(ord.data$housing)))) + 
        scale_color_manual(values=mgsub(c("wt", "ko"), 
                                        gt.cols, 
                                        levels(ord.data$host.gt))) + 
        facet_wrap(~ partition) + 
        theme(plot.margin=unit(c(0.01,0.01,0.01,0.01), "npc"),
              plot.background=element_blank(),
              strip.background=element_blank(), 
              strip.text=element_text(hjust=0)) + 
        guides(color=FALSE, shape=FALSE)
    partOrd.list[[i]] <- save.ord.plot
    
    ord.data$All <- rep("all", nrow(ord.data))
    beta.disp.data <- data.frame()
    nmds.coords <- ord.data[,c("NMDS1", "NMDS2")]
    
    sub.data <- cbind(nmds.coords, 
                      melt(ord.data, 
                           id.vars="partition", 
                           measure.var="All")[, c("variable", "value", "partition")])
    ordiplot(ord, display="sites")
    spider <- with(sub.data, ordispider(ord,
                                        display="sites",
                                        groups=value,
                                        spiders="centroid"))
    spider.df <- as.data.frame(spider[,])
    names(spider.df) <- c("cent.NMDS1", "cent.NMDS2")
    ord.coords <- sub.data[,c(1:2)]
    all.coords <- merge(ord.coords, spider.df, by="row.names")
    dist.to.centroid <- apply(all.coords, 1, function(x) euclid.dist(x))
    sub.data <- cbind(sub.data,
                      dist.to.centroid)
    beta.disp.data <- rbind(beta.disp.data, sub.data)
    
    factor <- "host.gt"
    sub.data <- cbind(nmds.coords, 
                      melt(ord.data, 
                           id.vars="partition",
                           measure.var=factor)[, c("variable", "value", "partition")])
    ordiplot(ord, display="sites")
    spider <- with(sub.data, ordispider(ord,
                                        display="sites",
                                        groups=value,
                                        spiders="centroid"))
    spider.df <- as.data.frame(spider[,])
    names(spider.df) <- c("cent.NMDS1", "cent.NMDS2")
    ord.coords <- sub.data[,c(1:2)]
    all.coords <- merge(ord.coords, spider.df, by="row.names")
    dist.to.centroid <- apply(all.coords, 1, function(x) euclid.dist(x))
    sub.data <- cbind(sub.data,
                      dist.to.centroid)
    beta.disp.data <- rbind(beta.disp.data, sub.data)
    betadisp.list[[i]] <- beta.disp.data
}
legend <- ord.plot +
    geom_point(aes(color=host.gt, shape=housing), size=2) + 
    scale_shape_manual(name="Housing", 
                       values=housing.shps,
                       labels=mgsub(c("isolated", "cohoused"), 
                                    c("Segregated", "Mixed"), 
                                    levels(ord.data$housing))) + 
    scale_color_manual(name="Host genotype", 
                       values=gt.cols, 
                       labels=c("Wild type", 
                                expression(italic(paste("rag1"))^-{}))) + 
    coord_cartesian(xlim = c(9998,9999)) + 
    theme(legend.position=c(0.5, 0.5), 
          legend.box="horizontal",
          legend.text.align=0, 
          plot.margin=unit(c(0.01,0.01,0.01,0.01), "npc"),
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.line=element_blank(),
          axis.ticks=element_blank(),
          panel.border=element_blank(),
          panel.grid=element_blank())

partitioned.ords <- arrangeGrob(partOrd.list[[1]], partOrd.list[[2]], partOrd.list[[3]], legend, 
                                ncol=2, as.table=FALSE)
ggsave(partitioned.ords,
       file="B_Plots/paritioned_ordinations.pdf",
       units="mm",
       width=86,
       height=86)

############

# PERMANOVA Stats on ordinations
load(partitioned_physeq_file, verbose=TRUE)

dist.obj <- phyloseq::distance(partitioned.physeq, method=dist.method, binary=dist.binary, na.rm=TRUE)
partitioned.data <- as(sample_data(partitioned.physeq), "data.frame")

adns <- adonis(dist.obj ~ partition * housing * host.gt,
               data=partitioned.data,
               permutations=9999)
save(adns, file="E_Saved_objects/partitioned_permanova_full")

# load("E_Saved_objects/partitioned_permanova_full", verbose=TRUE)
adns.stats <- summary(permustats(adns))
adns.ses <- merge(adns$aov.tab, adns.stats$z, by=0, sort=FALSE, all=TRUE)
adns.ses.tbl <- data.frame("Df"=adns.ses$Df,
                           "SumOfSqs" = adns.ses$SumsOfSqs,
                           "MeanSqs" = adns.ses$MeanSqs,
                           "F" = adns.ses$F.Model,
                           "R2"=adns.ses$R2, 
                           `Pr(>F)`=adns.ses$`Pr(>F)`,
                           "SES"=adns.ses$y,
                           row.names=mgsub(c("partition", "housing", "host.gt"),
                                           c("Partition", "Housing", "Host Genotype"),
                                           adns.ses$Row.names))
adns.ses.tbl
write.table(adns.ses.tbl, 
            file=paste("C_Text_Tables/allParts_ord_adns_",
                       mgsub(c(" ", "ø"), c("", "o"), dist.name),
                       ".txt",
                       sep=""), 
            sep="\t", 
            quote=FALSE)

# Custom PERMANOVA contrasts (ad hoc post hoc tests)
# Below are the dummy codes I generated for each comparison. (Not all are used in the final adonis())
host.gt <- partitioned.data$host.gt
housing <- partitioned.data$housing
partition <- partitioned.data$partition
contrasts(host.gt) <- c(-1, 1)
Host.gt <- model.matrix(~ host.gt)[, -1]
contrasts(housing) <- c(-1, 1)
Housing <- model.matrix(~ housing)[, -1]
contrasts(partition) <- cbind(c(0,1,0), c(0,0,1))
host.gt.in.neutral <- Host.gt * ifelse(partition == "neutral", 1, 0)
host.gt.in.overRep <- Host.gt * ifelse(partition == "overRep", 1, 0)
host.gt.in.underRep <- Host.gt * ifelse(partition == "underRep", 1, 0)
housing.in.neutral <- Housing * ifelse(partition == "neutral", 1, 0)
housing.in.overRep <- Housing * ifelse(partition == "overRep", 1, 0)
housing.in.underRep <- Housing * ifelse(partition == "underRep", 1, 0)
host.gt.in.neutral.iso <- Host.gt * ifelse(partition == "neutral" &
                                               housing == "isolated", 1, 0)
host.gt.in.neutral.coh <- Host.gt * ifelse(partition == "neutral" &
                                               housing == "cohoused", 1, 0)
host.gt.in.overRep.iso <- Host.gt * ifelse(partition == "overRep" &
                                               housing == "isolated", 1, 0)
host.gt.in.overRep.coh <- Host.gt * ifelse(partition == "overRep" &
                                               housing == "cohoused", 1, 0)
host.gt.in.underRep.iso <- Host.gt * ifelse(partition == "underRep" &
                                                housing == "isolated", 1, 0)
host.gt.in.underRep.coh <- Host.gt * ifelse(partition == "underRep" &
                                                housing == "cohoused", 1, 0)
contrasts(partition) <- cbind(c(1, -1, 0), c(1, 0, -1))
Partition.ortho <- model.matrix(~ partition)[, -1]
neutralVoverRep <- Partition.ortho[, 1]
neutralVunderRep <- Partition.ortho[, 2]

# This PERMANOVA takes almost 4 minutes to run on my computer
custom.adns <- adonis(dist.obj ~
                          host.gt.in.neutral.iso + host.gt.in.neutral.coh +
                          host.gt.in.overRep.iso + host.gt.in.overRep.coh +
                          host.gt.in.underRep.iso + host.gt.in.underRep.coh,
                      perm=9999)
save(custom.adns, file="E_Saved_objects/partitioned_permanova_custom_contrasts")

# load("E_Saved_objects/partitioned_permanova_custom_contrasts", verbose=TRUE)
custom.adns.stats <- summary(permustats(custom.adns))
custom.adns.ses <- merge(custom.adns$aov.tab, custom.adns.stats$z, by=0, sort=FALSE, all=TRUE)
custom.adns.ses.tbl <- data.frame("Df"=custom.adns.ses$Df, 
                                  "SumOfSqs" = custom.adns.ses$SumsOfSqs,
                                  "MeanSqs" = custom.adns.ses$MeanSqs,
                                  "F" = custom.adns.ses$F.Model,
                                  "R2"=custom.adns.ses$R2, 
                                  `Pr(>F)`=custom.adns.ses$`Pr(>F)`, 
                                  "SES"=custom.adns.ses$y,
                                  row.names=mgsub(c("neutral", 
                                                    "V", 
                                                    "overRep", 
                                                    "underRep", 
                                                    "housing", 
                                                    "host.gt", 
                                                    ".in.", 
                                                    ".iso", 
                                                    ".coh"), 
                                                  c("Neutral", 
                                                    " vs. ", 
                                                    "Over represented", 
                                                    "Under represented", 
                                                    "Housing", 
                                                    "Host Genotype", 
                                                    " in ", 
                                                    "-Segregated", 
                                                    "-Mixed"), 
                                                  custom.adns.ses$Row.names))
custom.adns.ses.tbl
write.table(custom.adns.ses.tbl, 
            file=paste("C_Text_Tables/allParts_ord_custAdns_", 
                       mgsub(c(" ", "ø"), c("", "o"), dist.name), 
                       ".txt", 
                       sep=""), 
            sep="\t", 
            quote=FALSE)
