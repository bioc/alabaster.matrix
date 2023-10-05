#' Array loading utilities
#'
#' Utilities for loading an array saved by \code{\link{stageObject}}. 
#'
#' @param info A named list of metadata for this array.
#' @param project Any argument accepted by the acquisition functions, see \code{?\link{acquireFile}}. 
#' By default, this should be a string containing the path to a staging directory.
#' @param names Logical scalar indicating whether the seed should be annotated with dimnames (if available).
#' @param path String containing the path to the file containing said array.
#' @param group String containing the name of the group with the dimnames.
#' @param ndim Integer scalar specifying the number of dimensions.
#'
#' @return \code{.createRawArraySeed} returns a seed that can be used in the \code{\link{DelayedArray}} constructor.
#'
#' \code{.extractArrayDimnames} returns a list of character vectors or \code{NULL}, containing the dimnames.
#'
#' @details
#' For \code{.createArraySeed}, the array should be one of:
#' \itemize{
#' \item \code{hdf5_dense_array}
#' \item \code{hdf5_sparse_matrix}
#' \item \code{hdf5_delayed_array}
#' \item \code{amalgamated_array}
#' }
#'
#' For delayed arrays, the file may contain a seed array with the \code{"custom alabaster local array"} type.
#' This should have a \code{path} dataset containing a relative path to another array in the same \code{project}, which is loaded and used as the seed for this delayed array.
#' Callers can overwrite this behavior by setting \code{"custom alabaster local array"} in the \code{knownArrays} from \pkg{chihaya} before calling \code{.createRawArraySeed}.
#' 
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
#' .createRawArraySeed(meta, project=dir)
#'
#' @export
#' @name createRawArraySeed
#' @importFrom HDF5Array HDF5ArraySeed H5SparseMatrixSeed
#' @importFrom rhdf5 h5read h5readAttributes
#' @importFrom alabaster.base loadObject acquireFile acquireMetadata
.createRawArraySeed <- function(info, project, names = TRUE) {
    path <- acquireFile(project, info$path)

    if ("hdf5_delayed_array" %in% names(info)) {
        copy <- chihaya::knownArrays()
        key <- "custom alabaster local array"
        if (!(key %in% names(copy))) {
            copy[[key]] <- function(file, name, contents) {
                location <- h5read(file, file.path(name, "path"))
                info <- acquireMetadata(project, location)
                loadObject(info, project)
            }
            old <- chihaya::knownArrays(copy)
            on.exit(chihaya::knownArrays(old))
        }

        group <- info$hdf5_delayed_array$group
        return(chihaya::loadDelayed(path, group))
    }

    if ("hdf5_dense_array" %in% names(info)) {
        dtype <- if (identical(info$array$type, "boolean")) "logical" else NA
        ds <- info$hdf5_dense_array$dataset
        out <- HDF5ArraySeed(filepath=path, name=ds, type=dtype)

        attrs <- h5readAttributes(path, ds)
        miss.place <- attrs[["missing-value-placeholder"]]
        if (!is.null(miss.place)) {
            out <- DelayedMask(out, placeholder=miss.place)
        }

        name.group <- if (names) info$hdf5_dense_array$dimnames else NULL
        return(.array_namer(out, path, name.group))
    }

    if ("hdf5_sparse_matrix" %in% names(info)) {
        group <- info$hdf5_sparse_matrix$group
        out <- H5SparseMatrixSeed(filepath=path, group=group)

        attrs <- h5readAttributes(path, paste0(group, "/data"))
        miss.place <- attrs[["missing-value-placeholder"]]
        if (!is.null(miss.place)) {
            out <- DelayedMask(out, placeholder=miss.place)
        }

        if (identical(info$array$type, "boolean")) {
            out <- (DelayedArray(out) != 0L)@seed
        }

        name.group <- if (names) info$hdf5_sparse_matrix$dimnames else NULL
        return(.array_namer(out, path, name.group))
    }

    if ("amalgamated_array" %in% names(info)) {
        return(.load_amalgamated_array(info, project))
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
