staging.options <- new.env()
staging.options$preserve.delayed <- FALSE
staging.options$recycle.hdf5 <- FALSE

#' Preserve delayed operations during staging
#'
#' Preserve delayed operations via \pkg{chihaya} when staging a \linkS4class{DelayedArray} with \code{\link{stageObject}}.
#'
#' @param preserve Whether to preserve delayed operations using the \pkg{chihaya} specification.
#'
#' @return Logical scalar indicating whether delayed operations are to be preserved by the DelayedArray method.
#' If \code{preserve} is supplied, it is used to set this scalar, and the \emph{previous} value of the scalar is invisibly returned.
#'
#' @details
#' By default, any DelayedArray in \code{\link{stageObject}} will be saved as a new dense array or sparse matrix.
#' However, if this option is enabled, DelayedArrays will instead be saved in the \pkg{chihaya} specification,
#' where the delayed operations are themselves stored in the HDF5 file (see \url{https://artifactdb.github.io/chihaya/} for details).
#'
#' The \pkg{chihaya} specification is more complicated to parse but can be helpful in reducing disk usage.
#' One simple example is to avoid sparsity-breaking or integer-to-float operations by storing their delayed representations in the file.
#' If the seed matrix is derived from some immutable reference location, advanced users can even store links to that location instead of duplicating the seed data.
#'
#' @author Aaron Lun
#'
#' @examples
#' preserveDelayedOperations()
#' old <- preserveDelayedOperations(TRUE)
#' preserveDelayedOperations()
#' preserveDelayedOperations(old)
#' 
#' @export
preserveDelayedOperations <- function(preserve) {
    prev <- staging.options$preserve.delayed
    if (missing(preserve)) {
        prev
    } else {
        staging.options$preserve.delayed <- preserve
        invisible(prev)
    }
}

#' Recycle existing HDF5 files
#'
#' Re-use existing files in HDF5-backed arrays rather than reserializing them in \code{\link{stageObject}}.
#'
#' @param recycle Whether to recycle existing files for HDF5-backed DelayedArrays.
#' 
#' @return Logical scalar indicating whether HDF5 files are to be reused.
#' If \code{recycle} is supplied, it is used to set this scalar, and the \emph{previous} value of the scalar is invisibly returned.
#'
#' @details
#' If this options is enabled, \code{stageObject} will attempt to link/copy existing files for any HDF5-backed DelayedArray instances -
#' most specifically, \linkS4class{HDF5Array} objects and \linkS4class{H5SparseMatrix} objects using the 10X format.
#' This avoids re-serialization of the data for faster staging.
#' It also allows advanced users to add their own customizations into the HDF5 file during staging, as long as they do not interfere with \code{\link{loadArray}}.
#'
#' By default, this option is disabled as the properties of the existing file are not known in the general case.
#' In particular, the file might contain other groups/datasets that are irrelevant, and use up extra disk space if copied; or confidential, and should not be stored in the staging directory.
#' Users should only enable this option if they have full control over the generation and contents of the backing HDF5 files.
#'
#' Also note that any dimnames on \code{x} will be ignored during recycling.
#'
#' @author Aaron Lun
#'
#' @examples
#' recycleHdf5Files()
#' old <- recycleHdf5Files(TRUE)
#' recycleHdf5Files()
#' recycleHdf5Files(old)
#' 
#' @export
recycleHdf5Files <- function(recycle) {
    prev <- staging.options$recycle.hdf5
    if (missing(recycle)) {
        prev
    } else {
        staging.options$recycle.hdf5 <- recycle
        invisible(prev)
    }
}

