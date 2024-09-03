#' Read a delayed array from disk
#'
#' Read a delayed high-dimensional array from its on-disk representation.
#' This is usually not directly called by users, but is instead called by dispatch in \code{\link{readObject}}.
#'
#' @param path String containing a path to a directory, itself created by the \code{\link{saveObject}} method for a delayed array.
#' @param metadata Named list of metadata for this object, see \code{\link{readObject}} for more details.
#' @param delayed_array.reload.args Named list of arguments to be passed to \code{\link{reloadDelayedObject}}.
#' @param ... Further arguments, ignored.
#' 
#' @return A multi-dimensional array-like object.
#'
#' @seealso
#' \code{"\link{saveObject,DelayedArray-method}"}, to create the directory and its contents.
#'
#' \code{\link{reloadDelayedObject}}, for the methods to reload each delayed operation.
#'
#' @author Aaron Lun
#'
#' @examples
#' arr <- array(rpois(10000, 10), c(50, 20, 10))
#' dimnames(arr) <- list(
#'    paste0("GENE_", seq_len(nrow(arr))),
#'    letters[1:20],
#'    NULL
#' )
#'
#' dir <- tempfile()
#' saveObject(arr, dir)
#' readObject(dir)
#' 
#' @export
readDelayedArray <- function(path, metadata, delayed_array.reload.args=list(), ...) {
    fpath <- file.path(path, "array.h5")
    fhandle <- H5Fopen(fpath, "H5F_ACC_RDONLY")
    on.exit(H5Fclose(fhandle))
    out <- do.call(reloadDelayedObject, c(list(fhandle, "delayed_array"), delayed_array.reload.args))
    ReloadedArray(path=path, seed=out)
}
