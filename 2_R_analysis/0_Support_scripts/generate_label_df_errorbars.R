require(car)
require(DTK)
require(multcompView)
require(ggplot2)
require(plyr)
source("0_Support_scripts/factorially_paste_lists.R")
source("0_Support_scripts/factorially_paste_rows.R")

generate_label_df_errorbars <- function(orig.data, 
                                        ggplot.obj, 
                                        aov.eq, 
                                        compare.panels=FALSE, 
                                        ss.type=c(1, 3))
{
    if(length(ss.type) > 1) {
        ss.type <- ss.type[1]
    }
    aov.list <- list()
    ggplot.data <- ggplot.obj$data
    ggplot.build <- ggplot_build(ggplot.obj)
    for (i in c(1:length(ggplot.build$data))) {
        if("linetype" %in% names(ggplot.build$data[[i]])) {
            ggplot.build.data <- as.data.frame(ggplot.build$data[[i]])
        }
    }
    num.panels <-  max(as.numeric(ggplot.build$panel$layout$PANEL), na.rm=TRUE)
    measure <- gsub(" ", "", strsplit(aov.eq, "~", fixed=TRUE)[[1]][1])
    aov.factors <- strsplit(gsub(" ", "", strsplit(aov.eq, "~", fixed=TRUE)[[1]][2]), "*", fixed=TRUE)[[1]]
    if(length(aov.factors) < 2) {
        factor.levels <- levels(eval(parse(text=paste("orig.data$", aov.factors, sep=""))))
    } else {
        factor.levels <- sapply(orig.data[,aov.factors], levels)
    }
    
    if(is.list(factor.levels)==TRUE) {
        factor.levels <- as.list(factor.levels)
        contrasts.order <- factorially_paste_lists(factor.levels)
    } else {
        factor.levels <- as.data.frame(factor.levels)
        contrasts.order <- factorially_paste_rows(factor.levels)
    }
    ult.contrasts <- as.data.frame(do.call(paste, c(ggplot.data[aov.factors], sep=":")))
    ult.contrasts[,1] <- factor(ult.contrasts[,1], levels=contrasts.order)
    ult.contrasts.name <- paste(aov.factors, collapse=":")
    names(ult.contrasts) <- ult.contrasts.name
    ggplot.data <- cbind(ggplot.data, ult.contrasts)
    if (compare.panels==FALSE) {
        if (num.panels > 1) {
            labels.df <- NULL
            facet.vars <- as.data.frame(ggplot.build$panel$layout[, !names(ggplot.build$panel$layout) %in% c("PANEL", 
                                                                                                             "ROW", 
                                                                                                             "COL", 
                                                                                                             "SCALE_X", 
                                                                                                             "SCALE_Y", 
                                                                                                             "AXIS_X", 
                                                                                                             "AXIS_Y")])
            facet.names <- names(ggplot.build$panel$layout)[!names(ggplot.build$panel$layout) %in% c("PANEL", 
                                                                                                     "ROW", 
                                                                                                     "COL", 
                                                                                                     "SCALE_X", 
                                                                                                     "SCALE_Y", 
                                                                                                     "AXIS_X", 
                                                                                                     "AXIS_Y")]
            names(facet.vars) <- facet.names
            panel <- 0
            for (i in c(1:length(facet.vars[,1]))) {
                # i <- 1
                # print(i)
                sub.ggplot.data <- ggplot.data
                panel <- panel + 1
                for (facet.name in facet.names) {
                    sub.ggplot.data <- subset(sub.ggplot.data, 
                                              eval(parse(text=paste(facet.name, 
                                                                    ' == "', 
                                                                    facet.vars[i, facet.name], 
                                                                    '"', 
                                                                    sep=""))))
                }
                if(ss.type == 1) {
                    # print("Tukey Test")
                    sub.ggplot.aov <- aov(eval(parse(text=aov.eq)), data=sub.ggplot.data)
                    aov.list[[panel]] <- sub.ggplot.aov
                    sub.ggplot.tukey <- TukeyHSD(sub.ggplot.aov)
                    final.level <- names(sub.ggplot.tukey[length(sub.ggplot.tukey)])
                    final.level.tick <- paste("`", final.level, "`", sep="")
                    Tukey.levels <- sub.ggplot.tukey[[final.level]][,4]
                    if(length(names(Tukey.levels))==0) {
                        names(Tukey.levels) <- row.names(sub.ggplot.tukey[[final.level]])
                    }
                    Tukey.labels <- multcompLetters(Tukey.levels, reversed=TRUE)['Letters']
                    Tukey.labels$Letters <- Tukey.labels$Letters[contrasts.order]
                    Tukey.labels$Letters <- Tukey.labels$Letters[!is.na(Tukey.labels$Letters)]
                    plot.labels <- names(Tukey.labels[['Letters']])
                } else if (ss.type == 3) {
                    # print("DTK Test")
                    sub.ggplot.aov <- Anova(lm(eval(parse(text=aov.eq)), sub.ggplot.data), type=ss.type)
                    aov.list[[panel]] <- sub.ggplot.aov
                    dtk.test <- with(sub.ggplot.data, 
                                     eval(parse(text=paste("DTK.test(", measure, ", interaction(", 
                                                           paste(aov.factors, collapse=","), 
                                                           "))", 
                                                           sep=""))))
                    dtk.results <- as.data.frame(dtk.test[[2]])
                    dtk.results <- within(dtk.results, {
                        p.value <- ifelse((is.na(`Lower CI`) | is.na(`Upper CI`)), 1, 
                                          ifelse((`Lower CI` < 0 & `Upper CI` < 0) | 
                                                     (`Lower CI` > 0 & `Upper CI` > 0), 0.01, 1))
                    })
                    Tukey.levels <- dtk.results$p.value
                    names(Tukey.levels) <- gsub(".", ":", row.names(dtk.results), fixed=TRUE)
                    Tukey.labels <- multcompLetters(Tukey.levels, reversed=TRUE)['Letters']
                    Tukey.labels$Letters <- Tukey.labels$Letters[contrasts.order]
                    plot.labels <- names(Tukey.labels[['Letters']])
                    final.level.tick <- paste(aov.factors, collapse=":")
                } else {
                    stop("ss.type must be 1 or 3")
                }
                
                max.range <- max(eval(parse(text=paste("sub.ggplot.data$", measure, sep=""))), na.rm=TRUE) - 
                    min(eval(parse(text=paste("sub.ggplot.data$", measure, sep=""))), na.rm=TRUE)
                buffer1 <- 10*(max.range^2)
                buffer2 <- 0.1*max.range
                if (buffer1 < buffer2) {
                    buffer <- buffer1
                } else {
                    buffer <- buffer2
                }
                boxplot.df <- ddply(sub.ggplot.data, 
                                    final.level.tick, 
                                    function (x) {
                                        max(eval(parse(text=paste("x$", measure, sep=""))), na.rm=TRUE) + 
                                            buffer
                                    })
                sub.ggplot.build.data <- subset(ggplot.build.data, PANEL == panel)
                x.coords <- sort(sub.ggplot.build.data$x)
                
                plot.levels <- data.frame(plot.labels, 
                                          sub.ggplot.data[c(1:length(x.coords)), facet.names], 
                                          x.coords, 
                                          labels = Tukey.labels[['Letters']], 
                                          stringsAsFactors = FALSE)
                names(plot.levels) <- c("contrast", facet.names, "x.coords", "labels")
                facet.labels.df <- merge(plot.levels, 
                                         boxplot.df, 
                                         by.x = "contrast", 
                                         by.y = 1, 
                                         sort = FALSE)
                labels.df <- rbind(labels.df, facet.labels.df)
            }
        } else {
            panel <- 1
            if(ss.type == 1) {
                # print("Tukey Test")
                ggplot.aov <- aov(eval(parse(text=aov.eq)), data=ggplot.data)
                aov.list[[panel]] <- ggplot.aov
                ggplot.tukey <- TukeyHSD(ggplot.aov)
                final.level <- names(ggplot.tukey[length(ggplot.tukey)])
                final.level.tick <- paste("`", final.level, "`", sep="")
                Tukey.levels <- ggplot.tukey[[final.level]][,4]
                if(length(names(Tukey.levels))==0) {
                    names(Tukey.levels) <- row.names(ggplot.tukey[[final.level]])
                }
                Tukey.labels <- multcompLetters(Tukey.levels, reversed=TRUE)['Letters']
                Tukey.labels$Letters <- Tukey.labels$Letters[contrasts.order]
                plot.labels <- names(Tukey.labels[['Letters']])
            } else if (ss.type == 3) {
                # print("DTK Test")
                ggplot.aov <- Anova(lm(eval(parse(text=aov.eq)), ggplot.data), type=ss.type)
                aov.list[[panel]] <- ggplot.aov
                dtk.test <- with(ggplot.data, 
                                 eval(parse(text=paste("DTK.test(", measure, ", interaction(", 
                                                       paste(aov.factors, collapse=","), 
                                                       "))", 
                                                       sep=""))))
                dtk.results <- as.data.frame(dtk.test[[2]])
                dtk.results <- within(dtk.results, {
                    p.value <- ifelse((is.na(`Lower CI`) | is.na(`Upper CI`)), 1, 
                                      ifelse((`Lower CI` < 0 & `Upper CI` < 0) | 
                                                 (`Lower CI` > 0 & `Upper CI` > 0), 0.01, 1))
                })
                Tukey.levels <- dtk.results$p.value
                names(Tukey.levels) <- gsub(".", ":", row.names(dtk.results), fixed=TRUE)
                Tukey.labels <- multcompLetters(Tukey.levels, reversed=TRUE, Letters = )['Letters']
                Tukey.labels$Letters <- Tukey.labels$Letters[contrasts.order]
                plot.labels <- names(Tukey.labels[['Letters']])
                final.level.tick <- paste(aov.factors, collapse=":")
            } else {
                stop("ss.type must be 1 or 3")
            }
            max.range <- max(eval(parse(text=paste("ggplot.data$", measure, sep=""))), na.rm=TRUE) - 
                min(eval(parse(text=paste("ggplot.data$", measure, sep=""))), na.rm=TRUE)
            buffer1 <- 10*(max.range^2)
            buffer2 <- 0.1*max.range
            if (buffer1 < buffer2) {
                buffer <- buffer1
            } else {
                buffer <- buffer2
            }
            boxplot.df <- ddply(ggplot.data, 
                                final.level.tick, 
                                function (x) {
                                    max(eval(parse(text=paste("x$", measure, sep=""))), na.rm=TRUE) + 
                                        buffer
                                })
            x.coords <- sort(ggplot.build.data$x)
            
            plot.levels <- data.frame(plot.labels, 
                                      x.coords, 
                                      labels = Tukey.labels[['Letters']], 
                                      stringsAsFactors = FALSE)
            labels.df <- merge(plot.levels, 
                               boxplot.df, 
                               by=1, 
                               sort = TRUE)
        }
    } else {
        panel <- 1
        facet.vars <- as.data.frame(ggplot.build$panel$layout[, !names(ggplot.build$panel$layout) %in% c("PANEL", 
                                                                                                         "ROW", 
                                                                                                         "COL", 
                                                                                                         "SCALE_X", 
                                                                                                         "SCALE_Y", 
                                                                                                         "AXIS_X", 
                                                                                                         "AXIS_Y")])
        facet.names <- names(ggplot.build$panel$layout)[!names(ggplot.build$panel$layout) %in% c("PANEL", 
                                                                                                 "ROW", 
                                                                                                 "COL", 
                                                                                                 "SCALE_X", 
                                                                                                 "SCALE_Y", 
                                                                                                 "AXIS_X", 
                                                                                                 "AXIS_Y")]
        names(facet.vars) <- facet.names
        if(ss.type == 1) {
            # print("Tukey Test")
            ggplot.aov <- aov(eval(parse(text=aov.eq)), data=ggplot.data)
            aov.list[[panel]] <- ggplot.aov
            ggplot.tukey <- TukeyHSD(ggplot.aov)
            final.level <- names(ggplot.tukey[length(ggplot.tukey)])
            final.level.tick <- paste("`", final.level, "`", sep="")
            Tukey.levels <- ggplot.tukey[[final.level]][,4]
            if(length(names(Tukey.levels))==0) {
                names(Tukey.levels) <- row.names(ggplot.tukey[[final.level]])
            }
            Tukey.labels <- multcompLetters(Tukey.levels, reversed=TRUE)['Letters']
            Tukey.labels$Letters <- Tukey.labels$Letters[contrasts.order]
            plot.labels <- names(Tukey.labels[['Letters']])
        } else if (ss.type == 3) {
            # print("DTK Test")
            ggplot.aov <- Anova(lm(eval(parse(text=aov.eq)), ggplot.data), type=ss.type)
            aov.list[[panel]] <- ggplot.aov
            dtk.test <- with(ggplot.data, 
                             eval(parse(text=paste("DTK.test(", measure, ", interaction(", 
                                                   paste(aov.factors, collapse=","), 
                                                   "))", 
                                                   sep=""))))
            dtk.results <- as.data.frame(dtk.test[[2]])
            dtk.results <- within(dtk.results, {
                p.value <- ifelse((is.na(`Lower CI`) | is.na(`Upper CI`)), 1, 
                                  ifelse((`Lower CI` < 0 & `Upper CI` < 0) | 
                                             (`Lower CI` > 0 & `Upper CI` > 0), 0.01, 1))
            })
            Tukey.levels <- dtk.results$p.value
            names(Tukey.levels) <- gsub(".", ":", row.names(dtk.results), fixed=TRUE)
            Tukey.labels <- multcompLetters(Tukey.levels, reversed=TRUE)['Letters']
            Tukey.labels$Letters <- Tukey.labels$Letters[contrasts.order]
            plot.labels <- names(Tukey.labels[['Letters']])
            final.level.tick <- paste(aov.factors, collapse=":")
        } else {
            stop("ss.type must be 1 or 3")
        }
        
        max.range <- max(eval(parse(text=paste("ggplot.data$", measure, sep=""))), na.rm=TRUE) - 
            min(eval(parse(text=paste("ggplot.data$", measure, sep=""))), na.rm=TRUE)
        buffer1 <- 10*(max.range^2)
        buffer2 <- 0.1*max.range
        if (buffer1 < buffer2) {
            buffer <- buffer1
        } else {
            buffer <- buffer2
        }
        boxplot.df <- ddply(ggplot.data, 
                            final.level.tick, 
                            function (x) {
                                max(eval(parse(text=paste("x$", measure, sep=""))), na.rm=TRUE) + 
                                    buffer
                            })
        x.coords.panels <- data.frame(x=ggplot.build.data$x, PANEL=ggplot.build.data$PANEL)
        sorted.x.coords.panels <- x.coords.panels[with(x.coords.panels, order(PANEL, x)),]
        x.coords <- sorted.x.coords.panels$x
        panel.names <- facet.vars[ggplot.build.data$PANEL, ,drop=FALSE]
        
        plot.levels <- data.frame(plot.labels, 
                                  panel.names,
                                  x.coords, 
                                  labels = Tukey.labels[['Letters']], 
                                  stringsAsFactors = FALSE)
        labels.df <- merge(plot.levels, 
                           boxplot.df, 
                           by=1, 
                           sort = TRUE)
    }
    return(list("labels.df"=labels.df, "aov.list"=aov.list))
}
