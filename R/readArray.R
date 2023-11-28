#' Read a dense array from disk
#'
#' Read a dense high-dimensional array from its on-disk representation.
#'
#' @param path String containing a path to a directory, itself created by the \code{\link{saveObject}} method for a dense array.
#' @param array.output.type String specifying the output type for this function.
#' This can be \code{"array"} or \code{"ReloadedArray"} (the default).
#' @param ... Further arguments, ignored.
#' 
#' @return A multi-dimensional array-like object.
#'
#' @details
#' By default, a file-backed \linkS4class{ReloadedArray} is returned to save memory and to preserve the provenance of the object.
#' Users can set \code{array.output.type="array"} to force \code{readArray} to load all data into memory and return an ordinary R array instead.
#'
#' @seealso
#' \code{"\link{saveObject,array-method}"}, to create the directory and its contents.
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
#' readArray(dir)
#' 
#' @export
#' @importFrom HDF5Array HDF5Array
#' @importFrom DelayedArray type<-
readArray <- function(path, array.output.type=NULL, ...) {
    fpath <- file.path(path, "array.h5")

    details <- local({
        fhandle <- H5Fopen(fpath)
        on.exit(H5Fclose(fhandle), add=TRUE, after=FALSE)
        ghandle <- H5Gopen(fhandle, "dense_array")
        on.exit(H5Gclose(ghandle), add=TRUE, after=FALSE)

        type <- h5_read_attribute(ghandle, "type")
        transposed <- h5_read_attribute(ghandle, "transposed", check=TRUE, default=0L)

        dhandle <- H5Dopen(ghandle, "data")
        on.exit(H5Dclose(dhandle), add=TRUE, after=FALSE)
        placeholder <- h5_read_attribute(dhandle, "missing-value-placeholder", check=TRUE, default=NULL)

        ndims <- H5Sget_simple_extent_dims(H5Dget_space(dhandle))$rank
        names <- load_names(ghandle, ndims)

        list(type=type, names=names, transposed=(transposed != 0L), placeholder=placeholder)
    })

    out <- HDF5Array(filepath=fpath, name="dense_array/data")
    if (type(out) == "raw") { # ... so that placeholders are correctly substituted.
        type(out) <- "integer"
    }

    if (!is.null(details$names)) {
        dimnames(out) <- rev(details$names)
    }
    if (!details$transposed) { # All R-based HDF5 bindings will automatically transpose, so we need to untranspose if `transposed=FALSE`.
        out <- t(out)
    }
    if (!is.null(details$placeholder)) {
        out <- DelayedMask(out, placeholder=details$placeholder)
        out <- DelayedArray(out)
    }
    type(out) <- from_array_type(details$type)

    if (is.null(array.output.type)) {
        array.output.type <- "ReloadedArray"
    } else {
        array.output.type <- match.arg(array.output.type, c("array", "ReloadedArray"))
    }
    if (array.output.type == "array") {
        return(as.array(out))
    }
    ReloadedArray(path=path, seed=out)
}

##############################
######### OLD STUFF ##########
##############################

#' @export
loadArray <- function(info, project) {
    seed <- .createRawArraySeed(info, project=project, names=TRUE)
    DelayedArray(seed)
}
