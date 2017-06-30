# neutral_theory_data.R

require(qdap)
require(doBy)
require(RColorBrewer)
require(plyr)
require(multcompView)
require(vegan)
require(picante)
require(phyloseq)
source("0_Support_scripts/Supplemental_Code_1.R")

neutral_theory_data <- function (physeq.object, 
                                 by.tank=FALSE,
                                 source.pool=FALSE) 
{
    curr <- c("gut", "water", "d09", "d75", "ko", "wt", "cohoused", "isolated")
    new <- c("Gut", "Water", "9 dpf", "75 dpf", "Knockout", "Wild type", "Mixed", "Segregated")
    genotypes <- as.character(unique(sample_data(physeq.object)$host.gt))
    house.types <- as.character(unique(sample_data(physeq.object)$housing))
    
    agg.plot.data <- data.frame(NULL)
    num.smpls <- NULL
    
    for (genotype in genotypes) {
        # genotype <- genotypes[1]
        for (house.type in house.types) {
            # house.type <- house.types[2]
            # print(paste(day, genotype, house.type, collapse="-"))
            smpl.data <- as(sample_data(physeq.object), "data.frame")
            pre.subset <- subset(smpl.data, host.gt==genotype & 
                                     housing==house.type)
            if(source.pool==TRUE) {
                pool_otu_file <- paste("E_Saved_objects/neutral_theory_data_source_pool_", 
                                       genotype, "_", 
                                       house.type, 
                                       "_otu_tbl.txt", 
                                       sep="")
                write.table(otu_table(subset_samples(physeq.object, host.gt==genotype & 
                                                         housing==house.type)), 
                            file=pool_otu_file, 
                            sep="\t")
                pool.otu.tbl <- as.matrix(read.table(pool_otu_file))
            }
            if (by.tank == TRUE) {
                tanks <- as.character(unique(pre.subset$tank))
                for (tank.name in tanks) {
                    # tank.name <- tanks[2]
                    data.subset <- subset(pre.subset, tank==tank.name)
                    keep.smpls <- row.names(data.subset)
                    curr.subset <- prune_samples(keep.smpls, physeq.object)
                    curr.subset <- prune_taxa(taxa_sums(curr.subset) > 0, curr.subset)
                    curr.nsmpls <- nsamples(curr.subset)
                    num.smpls <- c(num.smpls, curr.nsmpls)
                    otu_file <- paste("E_Saved_objects/otu_tbl_", 
                                      genotype, "_", 
                                      house.type, "_", 
                                      tank.name, ".txt", 
                                      sep="")
                    tax_file <- paste("E_Saved_objects/tax_tbl_", 
                                      genotype, "_", 
                                      house.type, "_", 
                                      tank.name, ".txt", 
                                      sep="")
                    write.table(otu_table(curr.subset), file=otu_file, sep="\t")
                    write.table(tax_table(curr.subset), file=tax_file, sep="\t")
                    otu.tbl <- as.matrix(read.table(otu_file))
                    tax.tbl <- as.matrix(read.table(tax_file))
                    
                    if(source.pool==TRUE) {
                        fit.tbl <- sncm.fit(otu.tbl, pool=pool.otu.tbl, taxon=tax.tbl, stats=FALSE)
                    } else {
                        fit.tbl <- sncm.fit(otu.tbl, taxon=tax.tbl, stats=FALSE)
                    }
                    good.name <- mgsub(curr, 
                                       new, 
                                       paste(genotype, 
                                             house.type, 
                                             tank.name, 
                                             sep="-"))
                    fit.tbl.df <- cbind(otu=row.names(fit.tbl), 
                                        name=rep(good.name, length(fit.tbl$p)), 
                                        fit.tbl)
                    row.names(fit.tbl.df) <- NULL
                    agg.plot.data <- rbind(agg.plot.data, fit.tbl.df)
                }
            } else {
                data.subset <- pre.subset
                keep.smpls <- row.names(data.subset)
                curr.subset <- prune_samples(keep.smpls, physeq.object)
                curr.subset <- prune_taxa(taxa_sums(curr.subset) > 0, curr.subset)
                curr.nsmpls <- nsamples(curr.subset)
                num.smpls <- c(num.smpls, curr.nsmpls)
                otu_file <- paste("E_Saved_objects/otu_tbl_", 
                                  genotype, "_", 
                                  house.type, ".txt", 
                                  sep="")
                tax_file <- paste("E_Saved_objects/tax_tbl_", 
                                  genotype, "_", 
                                  house.type, ".txt", 
                                  sep="")
                write.table(otu_table(curr.subset), file=otu_file, sep="\t")
                write.table(tax_table(curr.subset), file=tax_file, sep="\t")
                otu.tbl <- as.matrix(read.table(otu_file))
                tax.tbl <- as.matrix(read.table(tax_file))
                
                if(source.pool==TRUE) {
                    fit.tbl <- sncm.fit(otu.tbl, pool=pool.otu.tbl, taxon=tax.tbl, stats=FALSE)
                } else {
                    fit.tbl <- sncm.fit(otu.tbl, taxon=tax.tbl, stats=FALSE)
                }
                good.name <- mgsub(curr, new, paste(genotype, house.type, sep="-"))
                fit.tbl.df <- cbind(otu=row.names(fit.tbl), 
                                    name=rep(good.name, length(fit.tbl$p)), 
                                    fit.tbl)
                row.names(fit.tbl.df) <- NULL
                agg.plot.data <- rbind(agg.plot.data, fit.tbl.df)
            }
        }
    }
    # }
    overRep <- agg.plot.data[agg.plot.data$freq > agg.plot.data$pred.upr,]
    overRep <- overRep[!is.na(overRep$name),]
    overRep <- cbind(overRep, selection=rep("Over represented", length(overRep$otu)))
    underRep <- agg.plot.data[agg.plot.data$freq < agg.plot.data$pred.lwr,]
    underRep <- underRep[!is.na(underRep$name),]
    underRep <- cbind(underRep, selection=rep("Under represented", length(underRep$otu)))
    non.neutral <- rbind(overRep, underRep)
    neutral <- agg.plot.data[(agg.plot.data$freq <= agg.plot.data$pred.upr & 
                                  agg.plot.data$freq >= agg.plot.data$pred.lwr),]
    min.nsmpls <- min(num.smpls)
    return(list("all"=agg.plot.data, 
                "nonNeutral"=non.neutral, 
                "neutral"=neutral, 
                "min.smpls"=min.nsmpls))
}

