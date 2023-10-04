#' Delayed masking
#'
#' Delayed masking of missing values, based on replacement of placeholder values with \code{NA}.
#' This allows missingness to be encoded in frameworks without the same concept of \code{NA} as R.
#'
#' @param x An existing \pkg{DelayedArray} seed. 
#' @param placeholder Placeholder value to replace with \code{NA}.
#' This should be of the same type as \code{\link{type}(x)}.
#' @param force Whether to forcibly create a DelayedMask if \code{placeholder} is already \code{NA}.
#'
#' @return A DelayedMask object, to be wrapped in a \code{\link{DelayedArray}}.
#' If \code{force=FALSE} and \code{placeholder} is already \code{NA}, \code{x} is directly returned.
#'
#' @author Aaron Lun
#'
#' @aliases
#' DelayedMask-class
#' dim,DelayedMask-method
#' dimnames,DelayedMask-method
#' chunkdim,DelayedMask-method
#' path,DelayedMask-method
#' is_sparse,DelayedMask-method
#' extract_array,DelayedMask-method
#' extract_sparse_array,DelayedMask-method
#' OLD_extract_sparse_array,DelayedMask-method
#'
#' @examples
#' original <- DelayedArray(matrix(rpois(40, lambda=2), ncol=5))
#' original
#' masked <- DelayedMask(original, 0)
#' DelayedArray(masked)
#'
#' @export
#' @import methods
#' @importClassesFrom DelayedArray DelayedArray
DelayedMask <- function(x, placeholder, force=FALSE) {
    if (!force && is.na(placeholder)) {
        return(x)
    }

    if (is(x, "DelayedArray")) {
        x <- x@seed
    }
    new("DelayedMask", seed=x, placeholder=placeholder)
}

#' @export
#' @importClassesFrom DelayedArray DelayedUnaryOp
setClass("DelayedMask", contains="DelayedUnaryOp", slots=c(seed="ANY", placeholder="ANY"))

#' @export
setMethod("dim", "DelayedMask", function(x) callGeneric(x@seed))

#' @export
setMethod("dimnames", "DelayedMask", function(x) callGeneric(x@seed))

#' @export
#' @importFrom DelayedArray is_sparse
setMethod("is_sparse", "DelayedMask", function(x) {
    if (is.finite(x@placeholder) && x@placeholder == 0) {
        return(FALSE)
    } else {
        return(callGeneric(x@seed))
    }
})

#' @export
#' @importFrom DelayedArray extract_array
setMethod("extract_array", "DelayedMask", function(x, index) {
    ans <- callGeneric(x@seed, index)
    .replace_missing(ans, x@placeholder)
})

#' @export
#' @importFrom SparseArray extract_sparse_array
setMethod("extract_sparse_array", "DelayedMask", function(x, index) {
    ans <- callGeneric(x@seed, index)
    if (is(ans, "COO_SparseArray")) {
        ans@nzdata <- .replace_missing(ans@nzdata, x@placeholder)
    } else {
        ans@SVT <- .replace_missing_svt(ans@SVT, length(dim(x)) - 1L, x@placeholder)
    }
    ans
})

#' @export
#' @importFrom DelayedArray OLD_extract_sparse_array
setMethod("OLD_extract_sparse_array", "DelayedMask", function(x, index) {
    ans <- callGeneric(x@seed, index)
    ans@nzdata <- .replace_missing(ans@nzdata, x@placeholder)
    ans
})

.replace_missing <- function(vec, placeholder) {
    if (is.na(placeholder)) {
        if (is.nan(placeholder)) {
            vec[is.nan(vec)] <- NA
        }
    } else {
        # Using which() here to avoid problems with existing NAs.
        vec[which(vec == placeholder)] <- NA
    }
    vec
}

.replace_missing_svt <- function(tree, dim, placeholder) {
    if (dim == 1L) {
        for (i in seq_along(tree)) {
            current <- tree[[i]]
            if (!is.null(current)) {
                tree[[i]][[2]] <- .replace_missing(tree[[i]][[2]], placeholder)
            }
        }
    } else {
        for (i in seq_along(tree)) {
            current <- tree[[i]]
            if (!is.null(current)) {
                tree[[i]] <- .replace_missing_svt(current, dim - 1L, placeholder)
            }
        }
    }
    tree
}
