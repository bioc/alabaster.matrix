#' Read a sparse matrix from disk
#'
#' Read a sparse matrix from its on-disk representation.
#' This is usually not directly called by users, but is instead called by dispatch in \code{\link{readObject}}.
#'
#' @param path String containing a path to a directory, itself created by the \code{\link{saveObject}} method for a spars matrix.
#' @param metadata Named list of metadata for this object, see \code{\link{readObject}} for more details.
#' @param ... Further arguments, ignored.
#' 
#' @return A sparse \linkS4class{ReloadedMatrix} object.
#' 
#' @seealso
#' \code{"\link{saveObject,sparseMatrix-method}"}, to create the directory and its contents.
#'
#' @author Aaron Lun
#'
#' @examples
#' mat <- Matrix::rsparsematrix(100, 200, density=0.2)
#' rownames(mat) <- paste0("GENE_", seq_len(nrow(mat)))
#' dir <- tempfile()
#' saveObject(mat, dir)
#' readObject(dir)
#' 
#' @export
#' @importFrom HDF5Array H5SparseMatrixSeed
#' @importFrom DelayedArray type<-
readSparseMatrix <- function(path, metadata, ...) {
    fpath <- file.path(path, "matrix.h5")
    name <- "compressed_sparse_matrix"

    details <- local({
        fhandle <- H5Fopen(fpath, flags="H5F_ACC_RDONLY")
        on.exit(H5Fclose(fhandle), add=TRUE, after=FALSE)
        ghandle <- H5Gopen(fhandle, name)
        on.exit(H5Gclose(ghandle), add=TRUE, after=FALSE)

        type <- h5_read_attribute(ghandle, "type")
        layout <- h5_read_attribute(ghandle, "layout")
        shape <- h5_read_vector(ghandle, "shape")

        dhandle <- H5Dopen(ghandle, "data")
        on.exit(H5Dclose(dhandle), add=TRUE, after=FALSE)
        placeholder <- h5_read_attribute(dhandle, "missing-value-placeholder", check=TRUE, default=NULL)

        names <- load_names(ghandle, 2L)
        list(type=type, shape=shape, names=names, placeholder=placeholder, layout=layout)
    })

    out <- DelayedArray(H5SparseMatrixSeed(filepath=fpath, group=name, dim=details$shape, sparse.layout=tolower(details$layout)))
    if (type(out) == "raw") { # ... so that placeholders are correctly substituted.
        type(out) <- "integer"
    }

    if (!is.null(details$names)) {
        dimnames(out) <- details$names
    }
    if (!is.null(details$placeholder)) {
        out <- DelayedMask(out, placeholder=details$placeholder)
        out <- DelayedArray(out)
    }
    intended.type <- from_array_type(details$type)
    if (type(out) != intended.type) {
        type(out) <- intended.type
    }

    ReloadedArray(path, seed=out)
}
