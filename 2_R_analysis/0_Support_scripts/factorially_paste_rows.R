# factorially_paste_rows.R

factorially_paste_rows <- function(df) {
    total.combos <- nrow(df) ** ncol(df)
    combo.rows <- NULL
    for(column in c(1:ncol(df))) {
        rep.rowByCol <- NULL
        for(i in c(1:(nrow(df) ** (column - 1)))) {
            for(row in c(1:nrow(df))) {
                row.vec <- rep(as.character(df[row, column]), 
                               total.combos/(nrow(df) ** column))
                rep.rowByCol <- c(rep.rowByCol, row.vec)
            }
        }
        combo.rows <- cbind(combo.rows, rep.rowByCol)
    }
    output.vector <- apply(combo.rows, 1, paste, collapse=":")
    return(output.vector)
}