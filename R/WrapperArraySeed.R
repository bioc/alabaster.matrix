#' DelayedArray wrapper seed
#'
#' @description
#' Virtual class for a DelayedArray wrapper seed.
#' This automatically forwards DelayedArray generic operations onto an internal seed class.
#' Concrete subclasses are expected to attach more provenance-tracking information,
#' while the internal seed handles the heavy lifting of data extraction, e.g., \linkS4class{H5SparseMatrixSeed} or \linkS4class{HDF5ArraySeed} objects.
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
#' OLD_extract_sparse_array,WrapperArraySeed-method
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
#' OLD_extract_sparse_array(foo, list(1:10, 1:10))
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
#' @importFrom DelayedArray path
setMethod("path", "WrapperArraySeed", function(object, ...) callGeneric(object@seed))

#' @export
#' @importFrom DelayedArray is_sparse
setMethod("is_sparse", "WrapperArraySeed", function(x) callGeneric(x@seed))

#' @export
#' @importFrom DelayedArray extract_array
setMethod("extract_array", "WrapperArraySeed", function(x, index) callGeneric(x@seed, index))

#' @export
#' @importFrom DelayedArray OLD_extract_sparse_array
setMethod("OLD_extract_sparse_array", "WrapperArraySeed", function(x, index) callGeneric(x@seed, index))
