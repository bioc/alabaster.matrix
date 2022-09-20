#' DelayedArray wrapper seed
#'
#' @description
#' Virtual class for a DelayedArray wrapper seed.
#' This automatically forwards DelayedArray generic operations onto an internal seed class.
#' Concrete subclasses are expected to attach more provenance-tracking information,
#' while the internal seed handles the heavy lifting of data extraction, e.g., \linkS4class{H5SparseMatrixSeed} or \linkS4class{HDF5ArraySeed} objects.
#'
#' Subclass developers can also create methods for the \code{loadWrapperArray} generic.
#' This should accept two arguments:
#' \itemize{
#' \item \code{meta}, a list containing metadata for the array.
#' \item \code{project}, an object specifying the project of interest.
#' This is the sole argument used for S4 dispatch.
#' }
#' It should then return an instance of a WrapperArray subclass that retains some provenance about the resource from which it was generated.
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
#' @importFrom DelayedArray extract_sparse_array
setMethod("extract_sparse_array", "WrapperArraySeed", function(x, index) callGeneric(x@seed, index))


