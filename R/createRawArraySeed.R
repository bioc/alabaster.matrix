#' Array loading utilities
#'
#' Utilities for loading an array saved by \code{\link{stageObject}}. 
#'
#' @param info A named list of metadata for this array.
#' @param ndim Integer scalar specifying the number of dimensions.
#' @param path String containing the path to the file containing said array.
#' @param group String containing the name of the group with the dimnames.
#' @param names Logical scalar indicating whether the seed should be annotated with dimnames (if available).
#'
#' @return \code{.createRawArraySeed} returns a seed that can be used in the \code{\link{DelayedArray}} constructor.
#'
#' \code{.extractArrayDimnames} returns a list of character vectors or \code{NULL}, containing the dimnames.
#'
#' @details
#' For \code{.extractArrayDimnames}, \code{path} is expected to be a HDF5 file with a group specified by \code{group}.
#' Each child of this group is a string dataset named after a (0-indexed) dimension, containing the names for that dimension.
#'
#' @author Aaron Lun
#' @examples
#' # Staging an array as an example:
#' dir <- tempfile()
#' dir.create(dir)
#' mat <- array(rpois(10000, 10), c(50, 20, 10))
#' meta <- stageObject(mat, dir, "whee")
#'
#' # Loading it back as a DelayedArray seed:
#' .createRawArraySeed(meta, file.path(dir, meta$path))
#'
#' @export
#' @name createRawArraySeed
#' @importFrom HDF5Array HDF5ArraySeed H5SparseMatrixSeed
#' @importFrom rhdf5 h5readAttributes
.createRawArraySeed <- function(info, path, names = TRUE) {
    if ("hdf5_delayed_array" %in% names(info)) {
        group <- info$hdf5_delayed_array$group
        return(chihaya::loadDelayed(path, group))
    }

    if ("hdf5_dense_array" %in% names(info)) {
        dtype <- if (identical(info$array$type, "boolean")) "logical" else NA
        ds <- info$hdf5_dense_array$dataset
        out <- HDF5ArraySeed(filepath=path, name=ds, type=dtype)

        # Handling NA values for character arrays.
        if (type(out) == "character") {
            attrs <- h5readAttributes(path, ds)
            miss.place <- attrs[["missing-value-placeholder"]]
            if (!is.null(miss.place)) {
                out <- DelayedArray(out)
                out[out == miss.place] <- NA_character_
                out <- out@seed
            }
        }

        name.group <- if (names) info$hdf5_dense_array$dimnames else NULL
        return(.array_namer(out, path, name.group))
    }

    if ("hdf5_sparse_matrix" %in% names(info)) {
        group <- info$hdf5_sparse_matrix$group
        out <- H5SparseMatrixSeed(filepath=path, group=group)

        if (identical(info$array$type, "boolean")) {
            out <- (DelayedArray(out) != 0L)@seed
        }

        name.group <- if (names) info$hdf5_sparse_matrix$dimnames else NULL
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
    all.names <- which(info$group == paste0("/", group)) 
    for (i in all.names) {
        name <- info$name[i]
        dimnames[[as.integer(name) + 1L]] <- as.character(h5read(path, file.path("names", name)))
    }

    if (all(vapply(dimnames, is.null, TRUE))) {
        NULL
    } else {
        dimnames
    }
}
