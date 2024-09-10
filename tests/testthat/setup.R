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

setClass("SuperSeed", slots=c(dim="integer"))
setMethod("type", "SuperSeed", function(x) "integer")
setMethod("dim", "SuperSeed", function(x) x@dim)
setMethod("extract_array", "SuperSeed", function(x, index) {
    nr <- if (is.null(index[[1]])) x@dim[1] else length(index[[1]])
    nc <- if (is.null(index[[2]])) x@dim[2] else length(index[[2]])
    matrix(0L, nr, nc)
})
