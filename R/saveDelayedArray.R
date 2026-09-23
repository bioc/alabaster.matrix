#' Save DelayedArrays to disk
#'
#' Save \link{DelayedArray} objects to their on-disk representation.
#'
#' @param x A \linkS4class{DelayedArray} object.
#' @param path String containing a path to a directory in which to save \code{x}.
#' @param DelayedArray.dispatch.pristine Logical scalar indicating whether to call the \code{\link{saveObject}} methods of seeds of pristine arrays. 
#' @param DelayedArray.preserve.ops Logical scalar indicating whether delayed operations should be preserved on-disk.
#' @param DelayedArray.force.external Logical scalar indicating to save the seeds of \code{x} as external arrays.
#' This is passed directly to \code{storeDelayedObject} as the \code{save.external.array=} argument, see \code{?\link{storeDelayedObject}} for details.
#' @param DelayedArray.store.args More named arguments to pass to \code{\link{storeDelayedObject}}.
#' @param array.dedup.session,array.dedup.action Arguments controlling deduplication of \code{x}, see \code{?"\link{saveObject,array-method}"} for details.
#' If \code{x} is not a duplicate of an existing object, these arguments will be passed to further methods as described for \code{...}.
#' @param ... Further arguments passed to \code{storeDelayedObject} as \code{external.save.args}, if the delayed operations are to be preserved;
#' otherwise, they are passed to \code{\link{saveObject,array-method}} or \code{\link{saveObject,sparseMatrix-method}}.
#'
#' @details
#' Supplying \code{array.dedup.session=} by itself is only guaranteed to deduplicate \code{x} itself and may not deduplicate its seeds.
#' Users should combine this with \code{DelayedArray.force.external=TRUE} to force seeds to be saved via \code{saveObject},
#' which exposes the seeds to the deduplication machinery in their respective \code{saveObject} methods.
#' Check out \code{?"\link{storeDelayedObject}"} for more details.
#'
#' For finer control over which delayed operations are preserved with \code{DelayedArray.preserve.ops=TRUE},
#' we can construct \code{x} that contains \linkS4class{ToBeRealizedArray} seeds.
#' Any ToBeRealizedArray instance will force the realization of its wrapped delayed operations into the typical dense array/sparse matrix file representations.
#' This is useful if we only want to preserve some of the delayed operations,
#' e.g., if we want to realize a DelayedSubset from a larger dataset while preserving the delayed arithmetic on that subset.
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
setMethod("saveObject", "DelayedArray", function(
    x,
    path,
    DelayedArray.dispatch.pristine=TRUE,
    DelayedArray.preserve.ops=FALSE,
    DelayedArray.force.external=FALSE,
    DelayedArray.store.args=list(),
    array.dedup.session=NULL,
    array.dedup.action=NULL,
    ...)
{
    if (DelayedArray.dispatch.pristine && isPristine(x)) {
        s <- seed(x)
        out <- try_altSaveObject(
            s,
            path=path,
            DelayedArray.dispatch.pristine=DelayedArray.dispatch.pristine,
            DelayedArray.preserve.ops=DelayedArray.preserve.ops,
            DelayedArray.store.args=DelayedArray.store.args,
            array.dedup.session=array.dedup.session,
            array.dedup.action=array.dedup.action,
            ...
        )
        if (out) {
            return(invisible(NULL))
        }
    }

    if (!DelayedArray.preserve.ops) {
        if (is_sparse(x)) {
            .save_compressed_sparse_matrix(x, path, array.dedup.session=array.dedup.session, array.dedup.action=array.dedup.action, ...)
        } else {
            .save_array(x, path, array.dedup.session=array.dedup.session, array.dedup.action=array.dedup.action, ...)
        }
        return(invisible(NULL))
    }

    # This needs to be after the .save_array, .save_compressed_sparse_matrix
    # calls, as we don't want to add the object to the dedup session and then
    # have it get picked up by checkObjectInDedupSession within those calls.
    if (!is.null(array.dedup.session)) {
        dedup.path <- checkObjectInDedupSession(x, array.dedup.session)
        if (!is.null(dedup.path)) {
            cloneDirectory(dedup.path, path, array.dedup.action)
            return(invisible(NULL))
        }
        addObjectToDedupSession(x, array.dedup.session, path)
    }

    dir.create(path)
    saveObjectFile(path, "delayed_array", list(delayed_array=list(version="1.0")))

    fhandle <- H5Fcreate(file.path(path, "array.h5"))
    on.exit(H5Fclose(fhandle), add=TRUE, after=FALSE)

    if (!("external.save.args" %in% names(DelayedArray.store.args))) {
        # Technically, it seems that we should pass along the various
        # DelayedArray.* arguments for any call to altSaveObject to save
        # external seeds. However, there's no point, because any call to
        # saveObject,DelayedArray-method from storeDelayedObject will not
        # preserve delayed ops, otherwise we'd get an infinite recursion.
        DelayedArray.store.args$external.save.args <- list(..., array.dedup.session=array.dedup.session, array.dedup.action=array.dedup.action)
    }
    if (!("save.external.array" %in% names(DelayedArray.store.args))) {
        DelayedArray.store.args$save.external.array <- DelayedArray.force.external
    }
    do.call(altStoreDelayedObject, c(list(x@seed, handle=fhandle, name="delayed_array"), DelayedArray.store.args))

    ghandle <- H5Gopen(fhandle, "delayed_array")
    on.exit(H5Gclose(ghandle), add=TRUE, after=FALSE)
    h5_write_attribute(ghandle, "delayed_version", "1.1")

    invisible(NULL)
})
