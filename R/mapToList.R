#' Convert map from data.frame or DataFrame to list and vice versa
#'
#' The `mapToList` function provides a convenient way of reordering a
#' \code{data.frame} to a `list`. The `listToMap` function does the
#' opposite by taking a `list` and converting it to `DataFrame`.
#'
#' @param dfmap A \code{data.frame} or `DataFrame` object with
#' identifiers in the first column
#' @param assayCol A character vector of length one indicating the assay
#' names column
#' @return A `list` object of DataFrames for each assay
#' @example inst/scripts/listToMap-Ex.R
#' @export mapToList
mapToList <- function(dfmap, assayCol = "assay") {
    if (!S4Vectors::isSingleString(assayCol))
        stop("assay column name must be a single string")
    if (!assayCol %in% colnames(dfmap))
        stop("assay column not found in dataframe")
    assayColIndex <- which(names(dfmap) == assayCol)
    if (is.data.frame(dfmap))
        dfmap <- as(dfmap, "DataFrame")

    grp <- dfmap[[assayCol]]
    levs <- if (!length(grp)) levels(grp) else unique(grp)
    grp <- factor(grp, levels=levs, ordered = TRUE)

    S4Vectors::splitAsList(dfmap[, -assayColIndex], grp)
}
