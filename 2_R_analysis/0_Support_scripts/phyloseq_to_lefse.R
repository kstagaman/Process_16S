require(phyloseq)
options(scipen=999)

phyloseq_to_lefse <- function(phyloseq.object, 
                              row_names = c("housing", "host.gt", "smpl"), 
                              name = "phyloseq_object", 
                              write = TRUE, 
                              taxa.levels = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"), 
                              transpose_otus = FALSE) 
{
    smpl.data <- sample_data(phyloseq.object) # grab the mapping file from the phyloseq object
    smpl <- row.names(smpl.data)
    smpl.data <- cbind(smpl.data, smpl)
    t.smpl.data <- t(smpl.data)
    t.smpl.data <- as.data.frame(t.smpl.data) # transformed to put in proper structure
    
    t.smpl.data <- t.smpl.data[row_names, ] # only include rows of interest
    if (transpose_otus == TRUE) {
        otu.tbl <- t(as.matrix(otu_table(phyloseq.object))) # grab the otu table from the phyloseq object
    } else {
        otu.tbl <- otu_table(phyloseq.object) # grab the otu table from the phyloseq object
    }
    tax.tbl <- as(tax_table(phyloseq.object), "matrix") # grab the taxa table from the phyloseq object and coerce into a matrix
    tax.tbl <- tax.tbl[, taxa.levels]

# The following loop goes through the taxa table starting from the Domain level and gradually increasing to the lowest taxonomic level. For each level it appends the unique taxonomic levels to the uniq.lvls vector.
    uniq.lvls <- c()
    for(i in c(1:length(tax.tbl[1,]))) {
        lvls <- as.data.frame(do.call(paste, c(as.data.frame(tax.tbl[,1:i]), sep="|")))
        names(lvls) <- "tax.lvl"
        uniq.i <- as.character(unique(lvls$tax.lvl))
        uniq.lvls <- c(uniq.lvls, uniq.i)
    }
    tax.tbl.join <- as.data.frame(do.call(paste, c(as.data.frame(tax.tbl), sep="|")))
    row.names(tax.tbl.join) <- row.names(tax.tbl)
    names(tax.tbl.join) <- "tax.lvl"

# This loop goes through each sample (which are now column names for t.smpl.data), and calculates the relative abundance for each unique taxonomic levels (from above). These abundances only sum to 1 for *each taxonomic level*
    uniq.tax.lvl.abunds <- data.frame(row.names=uniq.lvls)
    for(smpl in names(t.smpl.data)) {
        abunds <- as.data.frame(otu.tbl[row.names(otu.tbl), smpl])
        total.abund <- sum(abunds[,1])
        smpl.tax.lvl.abunds <- cbind(abunds, tax.tbl.join)
        
        smpl.uniq.lvl.abunds <- data.frame()
        for(uniq.lvl in uniq.lvls) {
            uniq.sub <- subset(smpl.tax.lvl.abunds, grepl(uniq.lvl, smpl.tax.lvl.abunds$tax.lvl, fixed=TRUE))
            lvl.total <- as.factor(sum(uniq.sub[,1])/total.abund)
            uniq.lvl.smpl <- data.frame(row.names=uniq.lvl, "sample"=lvl.total)
            names(uniq.lvl.smpl) <- smpl
            smpl.uniq.lvl.abunds <- rbind(smpl.uniq.lvl.abunds, uniq.lvl.smpl)
        }

        uniq.tax.lvl.abunds <- cbind(uniq.tax.lvl.abunds, smpl.uniq.lvl.abunds)
    }
    
    final.data <- rbind(t.smpl.data, uniq.tax.lvl.abunds)

    if(write == TRUE) {
        file.name <- paste(gsub(".", "_", name, fixed=TRUE), "_lefse_data.txt", sep="")
        write.table(final.data, file=file.name, col.names=FALSE, sep="\t", quote=FALSE)
    } else {
        return(final.data)
    }
}
