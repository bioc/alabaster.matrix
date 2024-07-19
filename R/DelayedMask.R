#' Delayed masking
#'
#' Delayed masking of missing values, based on replacement of placeholder values with \code{NA}.
#' This allows missingness to be encoded in frameworks without the same concept of \code{NA} as R.
#'
#' @param x An existing \pkg{DelayedArray} seed. 
#' @param placeholder Placeholder value to replace with \code{NA}.
#' This should be of the same type as \code{\link{type}(x)}.
#'
#' @details
#' If \code{\link{is.na}(placeholder)} is true for double-precision \code{x}, masking is performed for all values of \code{x} where \code{is.na} is true.
#' This includes both NaNs and NAs; no attempt is made to distinguish between the NaN payloads.
#'
#' Currently, an error is raised for any integer \code{x} that produces non-missing values of -2^31 without a placeholder of \code{NA_integer_}.
#' This is because R cannot distinguish the integer -2^31 from an integer-type NA.
#'
#' @return A DelayedMask object, to be wrapped in a \code{\link{DelayedArray}}.
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
DelayedMask <- function(x, placeholder) {
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
#' @importFrom S4Arrays is_sparse
setMethod("is_sparse", "DelayedMask", function(x) {
    if (is.finite(x@placeholder) && x@placeholder == 0) {
        return(FALSE)
    } else {
        return(callGeneric(x@seed))
    }
})

#' @export
#' @importFrom S4Arrays extract_array
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
        nzidx <- nzwhich(ans)
        ans[nzidx] <- .replace_missing(ans[nzidx], x@placeholder)
    }
    ans
})

.replace_missing <- function(vec, placeholder) {
    if (is.na(placeholder)) {
        if (anyNA(vec)) {
            if (is.double(vec)) {
                vec[is.na(vec)] <- NA # convert NaNs as well.
            } else {
                # no-op, everything is already NA.
            }
        }
    } else {
        if (anyNA(vec)) {
            if (is.integer(vec)) {
                # We don't want to incorrectly propagate -2^31 as NAs, so we
                # just throw an error here. Better to fail than to silently
                # modify the meaning of the data.
                stop("integer arrays containing both -2^31 and NA are not yet supported")
            } else {
                vec[which(vec == placeholder)] <- NA # "which()" avoids problems with existing NaNs.
            }
        } else {
            vec[vec == placeholder] <- NA
        }
    }
    vec
}

