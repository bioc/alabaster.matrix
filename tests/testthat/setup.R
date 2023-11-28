expect_identical_without_names <- function(x, y) {
    if (!is.null(dimnames(x))) {
        if (identical(dimnames(x), vector("list", length(dim(x))))) {
            dimnames(x) <- NULL
        }
    }
    if (!is.null(dimnames(y))) {
        if (identical(dimnames(y), vector("list", length(dim(x))))) {
            dimnames(y) <- NULL
        }
    }
    expect_identical(x, y)
}
