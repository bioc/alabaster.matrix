#' DelayedArray wrapper seed
#'
#' @description
#' The WrapperArraySeed is, as the name suggests, a virtual class for a DelayedArray wrapper seed.
#' This forwards most of the DelayedArray generic operations onto an internal seed class,
#' typically a \linkS4class{H5SparseMatrixSeed} or \linkS4class{HDF5ArraySeed} objects from \code{\link{readSparseMatrix}} or \code{\link{readArray}}.
#' Similarly, the WrapperArray is a virtual DelayedArray class that contains a WrapperArraySeed.
#'
#' If an \pkg{alabaster} application operates on large arrays, developers may can consider defining concrete subclasses of the WrapperArraySeed (and WrapperArray).
#' These subclasses can store application-specific provenance-tracking information that persist throughout the lifetime of the array.
#' Such information is most useful for optimizing \code{\link{saveObject}} calls, which can instruct the application to link to the existing array rather than creating a new file.
#' Check out the \link{ReloadedArraySeed} class for an example of this approach.
#'
#' @aliases
#' loadWrapperArray
#' WrapperArraySeed-class
#' dim,WrapperArraySeed-method
#' dimnames,WrapperArraySeed-method
#' chunkdim,WrapperArraySeed-method
#' path,WrapperArraySeed-method
#' is_sparse,WrapperArraySeed-method
#' extract_array,WrapperArraySeed-method
#' extract_sparse_array,WrapperArraySeed-method
#' WrapperArray-class
#' coerce,WrapperArray,dgCMatrix-method
#' coerce,WrapperArraySeed,dgCMatrix-method
#'
#' @examples
#' # Mocking up a concrete wrapper array class, which contains an
#' # extra 'foo_id' slot to track the provenance of the data.
#' setClass("FooArraySeed", contains="WrapperArraySeed",
#'     slots=c(seed="ANY", foo_id="character"))
#'
#' y <- Matrix::rsparsematrix(1000, 100, 0.01)
#' foo <- new("FooArraySeed", seed=y, foo_id="FOO.0001")
#'
#' dim(foo)
#' is_sparse(foo)
#' extract_array(foo, list(1:10, 1:10))
#' extract_sparse_array(foo, list(1:10, 1:10))
#' 
#' @name WrapperArraySeed
NULL

#' @export
setMethod("dim", "WrapperArraySeed", function(x) callGeneric(x@seed))

#' @export
setMethod("dimnames", "WrapperArraySeed", function(x) callGeneric(x@seed))

#' @export
#' @importFrom DelayedArray chunkdim
setMethod("chunkdim", "WrapperArraySeed", function(x) callGeneric(x@seed))

#' @export
#' @importFrom BiocGenerics path
setMethod("path", "WrapperArraySeed", function(object, ...) callGeneric(object@seed))

#' @export
#' @importFrom S4Arrays is_sparse
setMethod("is_sparse", "WrapperArraySeed", function(x) callGeneric(x@seed))

#' @export
#' @importFrom S4Arrays extract_array
setMethod("extract_array", "WrapperArraySeed", function(x, index) callGeneric(x@seed, index))

#' @export
#' @importFrom SparseArray extract_sparse_array
setMethod("extract_sparse_array", "WrapperArraySeed", function(x, index) callGeneric(x@seed, index))

##############################

# We define the coercion methods here to give us the opportunity to
# short-circuit block processing if a more efficient method exists,
# e.g., direct reading of HDF5 sparse contents into a dgCMatrix.

.coerce_to_target_raw <- function(seed, target) {
    # Trying to coerce the naked seed first, but if it doesn't have a suitable
    # method, then we wrap it in a DelayedArray and try again.
    coercer <- selectMethod("coerce", list(class(seed), target), optional=TRUE)
    if (!is.null(coercer)) {
        out <- try(coercer(seed), silent=TRUE)
        # Some DelayedOp classes indicate that they do have a method for
        # coercion to target but it fails, so we need to catch that case.
        if (!is(out, "try-error")) {
            return(out)
        }
    }
    as(DelayedArray(seed), target)
}

.coerce_to_target1 <- function(wrapper, target) .coerce_to_target_raw(wrapper@seed, target)

#' @export
setAs("WrapperArraySeed", "dgCMatrix", function(from) .coerce_to_target1(from, "dgCMatrix"))

# The WrapperArray class only exists to define the coercion methods so that
# they don't have to be individually defined by each concrete subclass.

.coerce_to_target2 <- function(wrapper, target) .coerce_to_target_raw(wrapper@seed@seed, target)

#' @export
setAs("WrapperArray", "dgCMatrix", function(from) .coerce_to_target2(from, "dgCMatrix"))
