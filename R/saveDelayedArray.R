#' Save DelayedArrays to disk
#'
#' Save \link{DelayedArray} objects to their on-disk representation.
#'
#' @param x A \linkS4class{DelayedArray} object.
#' @param path String containing a path to a directory in which to save \code{x}.
#' @param delayedarray.preserve.ops Logical scalar indicating whether delayed operations should be preserved on-disk.
#' @param ... Further arguments, ignored.
#'
#' @return
#' \code{x} is saved to \code{path} and \code{NULL} is invisibly returned.
#'
#' @author Aaron Lun
#' @examples
#' mat <- Matrix::rsparsematrix(100, 200, density=0.2)
#' rownames(mat) <- paste0("GENE_", seq_len(nrow(mat)))
#' dmat <- DelayedArray::DelayedArray(mat) * 1
#'
#' dir <- tempfile()
#' saveObject(dmat, dir)
#' list.files(dir)
#'
#' @name saveDelayedArray
#' @aliases 
#' stageObject,DelayedArray-method
#' stageObject,DelayedMatrix-method
NULL

#' @export
#' @rdname saveDelayedArray
#' @importFrom DelayedArray isPristine seed
setMethod("saveObject", "DelayedArray", function(x, path, delayedarray.preserve.ops=FALSE, ...) {
    if (isPristine(x)) {
        s <- seed(x)
        fun <- selectMethod("saveObject", class(s), optional=TRUE)
        if (!is.null(fun)) {
            return(fun(s, path, ...))
        }
    }

    if (!delayedarray.preserve.ops) {
        if (is_sparse(x)) {
            .save_compressed_sparse_matrix(x, path, ...)
        } else {
            .save_array(x, path, ...)
        }
    } else {
        stop("preservation of delayed operations is not currently supported")
    }

    invisible(NULL)
})
