#generate_distance_data
require(reshape2)

generate_distance_data <- function(dist.obj, 
                                   smpl.data, 
                                   measures.of.interest, 
                                   id.factors=NULL, 
                                   verbose=FALSE) {
  dist.obj <- as.matrix(dist.obj)
  dist.obj <- melt(dist.obj)
  names(dist.obj) <- c("smpl1", "smpl2", "distance")
  dist.obj <- subset(dist.obj, smpl1 != smpl2)
  dist.obj.data <- dist.obj
  
  measures.smpl1s <- smpl.data[as.character(dist.obj$smpl1), 
                               measures.of.interest[sapply(smpl.data[, measures.of.interest], is.factor)]]
  measures.smpl1s <- as.vector(do.call(paste, c(measures.smpl1s, sep="_")))
  measures.smpl2s <- smpl.data[as.character(dist.obj$smpl2), 
                               measures.of.interest[sapply(smpl.data[, measures.of.interest], is.factor)]]
  measures.smpl2s <- as.vector(do.call(paste, c(measures.smpl2s, sep="_")))
  measures.smpl.pairs <- data.frame(measures.smpl1s, measures.smpl2s)
  sort.measures.smpl.pairs <- apply(measures.smpl.pairs, 1, sort)
  comp.df <- data.frame(full.comp=apply(sort.measures.smpl.pairs, 2, paste, collapse="<>"))
  comp.df$full.comp <- factor(comp.df$full.comp)
  dist.obj.data <- cbind(dist.obj.data, comp.df)

  for (measure in measures.of.interest) {
    # measure <- measures.of.interest[1]
    
    if (verbose == TRUE) {print(paste(measure, "...", sep=""), quote=FALSE)}
    
    if(is.factor(eval(parse(text=paste("smpl.data$", measure, sep=""))))==TRUE) {
      comparison <- paste(measure, ".comp", sep="")
      measure.smpl1s <- smpl.data[as.character(dist.obj$smpl1), measure]
      measure.smpl2s <- smpl.data[as.character(dist.obj$smpl2), measure]
      measure.smpl.pairs <- data.frame(measure.smpl1s, measure.smpl2s)
      sort.measure.smpl.pairs <- apply(measure.smpl.pairs, 1, sort)
      comp.df <- data.frame(comp=apply(sort.measure.smpl.pairs, 2, paste, collapse="<>"))
      comp.df$comp <- factor(comp.df$comp)
      names(comp.df) <- comparison
      dist.obj.data <- cbind(dist.obj.data, comp.df)
    } else {
      comparison <- paste(measure, ".diff", sep="")
      measure.smpl1s <- smpl.data[as.character(dist.obj$smpl1), measure]
      measure.smpl2s <- smpl.data[as.character(dist.obj$smpl2), measure]
      comp.df <- data.frame("comp"=abs(measure.smpl1s - measure.smpl2s))
      names(comp.df) <- comparison
      dist.obj.data <- cbind(dist.obj.data, comp.df)
    }
    if (verbose == TRUE) {print("done", quote=FALSE)}
  }
  for (id.factor in id.factors) {
      # id.factor <- id.factors[1]
      
      if (verbose == TRUE) {print(paste(id.factor, "...", sep=""), quote=FALSE)}
      
      if(is.factor(eval(parse(text=paste("smpl.data$", id.factor, sep=""))))==TRUE) {
          comparison <- paste(id.factor, ".comp", sep="")
          id.factor.smpl1s <- smpl.data[as.character(dist.obj$smpl1), id.factor]
          id.factor.smpl2s <- smpl.data[as.character(dist.obj$smpl2), id.factor]
          id.factor.smpl.pairs <- data.frame(id.factor.smpl1s, id.factor.smpl2s)
          sort.id.factor.smpl.pairs <- apply(id.factor.smpl.pairs, 1, sort)
          comp.df <- data.frame(comp=apply(sort.id.factor.smpl.pairs, 2, paste, collapse="<>"))
          comp.df$comp <- factor(comp.df$comp)
          names(comp.df) <- comparison
          dist.obj.data <- cbind(dist.obj.data, comp.df)
      } else {
          comparison <- paste(id.factor, ".diff", sep="")
          id.factor.smpl1s <- smpl.data[as.character(dist.obj$smpl1), id.factor]
          id.factor.smpl2s <- smpl.data[as.character(dist.obj$smpl2), id.factor]
          comp.df <- data.frame("comp"=abs(id.factor.smpl1s - id.factor.smpl2s))
          names(comp.df) <- comparison
          dist.obj.data <- cbind(dist.obj.data, comp.df)
      }
      if (verbose == TRUE) {print("done", quote=FALSE)}
  }
  return(dist.obj.data)
}

