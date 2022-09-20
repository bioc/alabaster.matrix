#' Array loading utilities
#'
#' Utilities for loading an array saved by \code{\link{stageObject}}. 
#'
#' @param info A named list of metadata for this array.
#' @param ndim Integer scalar specifying the number of dimensions.
#' @param path String containing the path to the file containing said array.
#' @param names Logical scalar indicating whether the seed should be annotated with dimnames (if available).
#'
#' @return \code{.createRawArraySeed} returns a seed that can be used in the \code{\link{DelayedArray}} constructor.
#' For matrices, this will call \code{\link{.createRawMatrixSeed}}.
#'
#' \code{.extractArrayDimnames} returns a list of character vectors or \code{NULL}, containing the dimnames.
#'
#' @details
#' For \code{.extractArrayDimnames}, \code{path} is expected to be a HDF5 file with a \code{names} group.
#' Each child of this group is a string dataset named after a (0-indexed) dimension, containing the names for that dimension.
#'
#' @author Aaron Lun
#'
#' @export
#' @name createRawArraySeed
#' @importFrom HDF5Array HDF5ArraySeed H5SparseMatrixSeed
.createRawArraySeed <- function(info, path, names = TRUE) {
    if ("hdf5_delayed_array" %in% names(info)) {
        group <- info$hdf5_delayed_array$group
        return(chihaya::loadDelayed(path, group))
    }

    if ("hdf5_dense_array" %in% names(info)) {
        group <- info$hdf5_dense_array$group
        name.group <- if (names) info$hdf5_dense_array$dimnames else NULL
        out <- HDF5ArraySeed(filepath=path, name=group)
        return(.array_namer(out, path, name.group))
    }

    if ("hdf5_sparse_matrix" %in% names(info)) {
        group <- info$hdf5_sparse_array$group
        name.group <- if (names) info$hdf5_sparse_matrix$dimnames else NULL
        out <- H5SparseMatrixSeed(filepath=path, group=group)
        return(.array_namer(out, path, name.group))
    }

    stop("unsupported type '", info$`_extra`$type, "'")
}

#' @importFrom DelayedArray DelayedArray
.array_namer <- function(seed, path, names.group) {
    if (!is.null(names.group)) {
        dnames <- .extractArrayDimnames(path, names.group, length(dim(seed)))
        if (!is.null(dnames)) {
            mat <- DelayedArray(seed)
            dimnames(mat) <- dnames
            seed <- mat@seed
        }
    }
    seed
}

#' @export
#' @rdname createRawArraySeed
#' @importFrom rhdf5 h5read h5ls
.extractArrayDimnames <- function(path, group, ndim) {
    info <- h5ls(path)

    dimnames <- vector("list", ndim)
    all.names <- which(info$group == paste0("/", group, "/", "names"))
    for (i in all.names) {
        name <- info$name[i]
        dimnames[[as.integer(name)]] <- as.character(h5read(path, file.path("names", name)))
    }

    if (all(vapply(dimnames, is.null, TRUE))) {
        NULL
    } else {
        dimnames
    }
}
