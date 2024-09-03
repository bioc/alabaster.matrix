#' Save DelayedArrays to disk
#'
#' Save \link{DelayedArray} objects to their on-disk representation.
#'
#' @param x A \linkS4class{DelayedArray} object.
#' @param path String containing a path to a directory in which to save \code{x}.
#' @param DelayedArray.dispatch.pristine Logical scalar indicating whether to call the \code{\link{saveObject}} methods of seeds of pristine arrays. 
#' @param DelayedArray.preserve.ops Logical scalar indicating whether delayed operations should be preserved on-disk.
#' @param DelayedArray.store.args Named arguments to pass to \code{\link{storeDelayedObject}}.
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
#' saveObject(dmat, dir, delayed.preserve.ops=TRUE)
#' list.files(dir)
#'
#' @seealso
#' \code{\link{storeDelayedObject}}, for the methods to save each delayed operation.
#'
#' @name saveDelayedArray
#' @aliases 
#' stageObject,DelayedArray-method
#' stageObject,DelayedMatrix-method
NULL

#' @export
#' @rdname saveDelayedArray
#' @importFrom DelayedArray isPristine seed
setMethod("saveObject", "DelayedArray", function(x, path, DelayedArray.dispatch.pristine=TRUE, DelayedArray.preserve.ops=FALSE, DelayedArray.store.args=list(), ...) {
    if (DelayedArray.dispatch.pristine && isPristine(x)) {
        s <- seed(x)
        fun <- selectMethod("saveObject", class(s), optional=TRUE)
        if (!is.null(fun)) {
            return(fun(s, path, ...))
        }
    }

    if (!DelayedArray.preserve.ops) {
        if (is_sparse(x)) {
            .save_compressed_sparse_matrix(x, path, ...)
        } else {
            .save_array(x, path, ...)
        }
    } else {
        dir.create(path)
        saveObjectFile(path, "delayed_array", list(delayed_array=list(version="1.0")))

        fhandle <- H5Fcreate(file.path(path, "array.h5"))
        on.exit(H5Fclose(fhandle), add=TRUE, after=FALSE)

        if (!("external.save.args" %in% names(DelayedArray.store.args))) {
            DelayedArray.store.args$external.save.args <- list(...)
        }
        do.call(storeDelayedObject, c(list(x@seed, handle=fhandle, name="delayed_array"), DelayedArray.store.args))

        ghandle <- H5Gopen(fhandle, "delayed_array")
        on.exit(H5Gclose(ghandle), add=TRUE, after=FALSE)
        h5_write_attribute(ghandle, "delayed_version", "1.1")
    }

    invisible(NULL)
})
