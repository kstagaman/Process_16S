# factorially_paste_lists.R

factorially_paste_lists <- function(lst) {
    list.lens <- NULL
    for(i in c(1:length(lst))) {
        len.i <- length(lst[[i]])
        list.lens <- c(list.lens, len.i)
    }
    rep.factors <- c(1, list.lens)
    total.combos <- prod(list.lens)
    combo.rows <- NULL
    columns <- NULL
    for(column in c(1:length(lst))) {
        columns <- c(columns, column)
        rep.rowByCol <- NULL
        for(i in c(1:prod(rep.factors[1:column]))) {
            for(row in c(1:length(lst[[column]]))) {
                row.vec <- rep(as.character(lst[[column]][row]), 
                               total.combos/prod(list.lens[1:column]))
                rep.rowByCol <- c(rep.rowByCol, row.vec)
            }
        }
        combo.rows <- cbind(combo.rows, rep.rowByCol)
    }
    output.vector <- apply(combo.rows, 1, paste, collapse=":")
    return(output.vector)
}