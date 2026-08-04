#' Save HDF5-backed arrays to disk
#'
#' Save \link[HDF5Array]{HDF5Array} ond \link[HDF5Array]{H5SparseMatrix} bjects to an on-disk representation.
#'
#' @param x A \linkS4class{DelayedArray} object.
#' @param path String containing a path to a directory in which to save \code{x}.
#' @param array.dedup.session,array.dedup.action Arguments controlling deduplication of \code{x}, see \code{?"\link{saveObject,array-method}"} for details.
#' @param ... Further arguments, currently ignored.
#'
#' @details
#' These methods are provided solely to avoid \code{\link{saveObject}} dispatching to \code{\link{saveObject,DelayedArray-method}} upon encountering a HDF5-backed array.
#' If this dispatch occurs and \code{DelayedArray.preserve.ops=TRUE}, \code{x} will be saved as a DelayedArray on disk.
#' This is not incorrect per se but adds unnecessary files to the on-disk representation.
#' The methods listed here will avoid creating such files and save a bit of disk space.
#'
#' @return
#' \code{x} is saved to \code{path} and \code{NULL} is invisibly returned.
#'
#' @author Aaron Lun
#' @examples
#' library(HDF5Array)
#' mat <- as(array(rpois(1000, 10), c(50, 20)), "HDF5Array")
#' dir <- tempfile()
#' saveObject(mat, dir)
#' list.files(dir)
#'
#' @name saveHDF5Array
NULL

#' @export
#' @importClassesFrom HDF5Array HDF5Array
#' @rdname saveHDF5Array
setMethod("saveObject", "HDF5Array", function(x, path, array.dedup.session=NULL, array.dedup.action=NULL, ...) {
    .save_array(x, path, array.dedup.session=array.dedup.session, array.dedup.action=array.dedup.action, ...)
})

#' @export
#' @importClassesFrom HDF5Array H5SparseMatrix
#' @rdname saveHDF5Array
setMethod("saveObject", "H5SparseMatrix", function(x, path, array.dedup.session=NULL, array.dedup.action=NULL, ...) {
    .save_compressed_sparse_matrix(x, path, array.dedup.session=array.dedup.session, array.dedup.action=array.dedup.action, ...)
})
