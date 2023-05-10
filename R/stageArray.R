#' Stage a multi-dimensional array for upload
#'
#' Stage a high-dimensional array in preparation for upload to DataSetDB.
#'
#' @param x An array, almost always integer or numeric, though logical and character matrices are also supported.
#' Alternatively, a \linkS4class{DelayedArray} or any instance of a \linkS4class{Matrix} class.
#' @param dir String containing the path to the staging directory.
#' @param path String containing the relative path to a subdirectory inside the staging directory, in which \code{x} is to be saved.
#' @param child Logical scalar indicating whether \code{x} is a child of a larger object.
#' @param preserve Whether to preserve delayed operations using the \pkg{chihaya} specification.
#' @param recycle Whether to recycle existing files for HDF5-backed DelayedArrays.
#' 
#' @return
#' For the \code{stageObject} methods, the array is saved into a single file at \code{file.path(dir, path)}, possibly after appending an arbitrary file extension. 
#' A named list is returned, containing at least:
#' \itemize{
#' \item \code{$schema}, a string specifying the schema to use to validate the metadata.
#' \item \code{path}, a string containing the path to the file inside the subdirectory, containing the assay contents.
#' \item \code{is_child}, a logical scalar equal to the input \code{child}.
#' }
#'
#' For \code{preserveDelayedOperations}, a logical scalar is returned indicating whether delayed operations are to be preserved by the DelayedArray method.
#' If \code{preserve} is supplied, it is used to set this scalar, and the \emph{previous} value of the scalar is returned.
#'
#' @details
#' The default behavior is to save the array as a dense matrix in a HDF5 file using methods from the \pkg{HDF5Array} package.
#' Other representations may have more appropriate formats, which are supported by simply writing new methods for this generic.
#' Note that specialized methods will usually require new schemas to validate any new metadata fields.
#'
#' If \code{x} itself is a child of a larger object, we suggest using the output \code{path} when referencing \code{x} from within the larger object's metadata.
#' This is because \code{stageObject} methods may add more path components, file extensions, etc. to the input \code{path} when saving the object.
#' As a result, the output \code{path} may not be the same as the input \code{path}.
#'
#' By default, \code{preserveDelayedOperations()} is \code{FALSE} so any DelayedArray \code{x} will be saved as a dense HDF5 dataset.
#' If \code{preserveDelayedOperations()} is \code{TRUE}, DelayedArrays will instead be saved in the \pkg{chihaya} specification,
#' where the delayed operations are themselves stored in the HDF5 file (see \url{https://ltla.github.io/chihaya} for details).
#'
#' If \code{recycleHdf5Files()} is \code{TRUE}, \code{stageObject} will attempt to link/copy existing files for any HDF5-backed DelayedArray instances 
#' (most specifically, \linkS4class{HDF5Array} objects and \linkS4class{H5SparseMatrix} objects using the 10X format).
#' This avoids re-serialization of the data for faster staging, at the risk of accidentally including other data in the existing file that was not meant to be staged.
#' Users should only enable this option if it is known that the existing file contains nothing but the matrix data.
#' Also note that any dimnames on \code{x} will be ignored during recycling.
#'
#' @author Aaron Lun
#' @examples
#' dir <- tempfile()
#' dir.create(dir)
#'
#' mat <- array(rpois(10000, 10), c(50, 20, 10))
#' dimnames(mat) <- list(
#'    paste0("GENE_", seq_len(nrow(mat))),
#'    letters[1:20],
#'    NULL
#' )
#'
#' path <- "whee"
#' stageObject(mat, dir, path)
#'
#' list.files(dir)
#'
#' @name stageArray
#' @importFrom alabaster.base stageObject
NULL

#' @importFrom DelayedArray is_sparse
#' @importFrom rhdf5 h5createFile 
#' @importFrom HDF5Array writeHDF5Array
#' @importFrom alabaster.base .chooseMissingStringPlaceholder .addMissingStringPlaceholderAttribute
.stage_array <- function(x, dir, path, child=FALSE) {
    dir.create(file.path(dir, path), showWarnings=FALSE)
    xpath <- file.path(path, "array.h5")
    ofile <- file.path(dir, xpath)

    missing.placeholder <- NULL
    if (type(x) == "character" && anyNA(x)) {
        missing.placeholder <- .chooseMissingStringPlaceholder(x)
        x[is.na(x)] <- missing.placeholder
    }

    h5createFile(ofile)
    writeHDF5Array(x, filepath=ofile, name="data")
    nm <- .name_saver(x, ofile)

    if (!is.null(missing.placeholder)) {
        .addMissingStringPlaceholderAttribute(ofile, "data", missing.placeholder)
    }

    list(
        `$schema` = "hdf5_dense_array/v1.json",
        path = xpath,
        is_child = child,
        `array` = .grab_array_type(x),
        hdf5_dense_array = list(
            dataset = "data",
            dimnames = nm
        )
    )
}

.grab_array_type <- function(x) {
    list(
        dimensions=I(dim(x)),
        type=array_type(x)
    )
}

#' @importFrom rhdf5 h5createGroup h5write
.name_saver <- function(x, path, group = "names") {
    # Chucking in some names. 
    if (!is.null(dimnames(x))) {
        h5createGroup(path, group)
        d <- dimnames(x)
        for (i in seq_along(d)) {
            current <- d[[i]]
            if (!is.null(current)) {
                h5write(current, file=path, name=paste0(group, "/", i - 1L)) # 0-based
            }
        }
        return(group)
    } else {
        return(NULL)
    }
}

#' @export
#' @rdname stageArray
setMethod("stageObject", "array", function(x, dir, path, child=FALSE) .stage_array(x, dir, path, child=child))

#' @importFrom alabaster.base .stageObject
.stage_delayed <- function(x, dir, path, child, fallback) {
    out <- .check_for_hdf5(x, dir, path, child = child)
    if (!is.null(out)) {
        return(out)
    }

    if (!preserveDelayedOperations()) {
        return(fallback(x, dir, path, child=child))
    }

    dir.create(file.path(dir, path), showWarnings=FALSE)
    xpath <- file.path(path, "delayed.h5")
    chihaya::saveDelayed(x, file.path(dir, xpath), "data")

    list(
        `$schema`="hdf5_delayed_array/v1.json",
        path=xpath,
        is_child=child,
        `array` = .grab_array_type(x),
        hdf5_delayed_array= list(
            group="data"
        )
    )
}

#' @export
#' @rdname stageArray
setMethod("stageObject", "DelayedArray", function(x, dir, path, child=FALSE) .stage_delayed(x, dir, path, child = child, fallback = .stage_array))

staging.options <- new.env()
staging.options$preserve.delayed <- FALSE
staging.options$recycle.hdf5 <- FALSE

#' @export
#' @rdname stageArray
preserveDelayedOperations <- function(preserve) {
    prev <- staging.options$preserve.delayed
    if (missing(preserve)) {
        prev
    } else {
        staging.options$preserve.delayed <- preserve
        invisible(prev)
    }
}

#' @importFrom rhdf5 h5createFile h5createGroup
#' @importFrom HDF5Array writeHDF5Array
.stage_sparse_matrix <- function(x, dir, path, child=FALSE) {
    dir.create(file.path(dir, path), showWarnings=FALSE)
    xpath <- file.path(path, "matrix.h5")
    ofile <- file.path(dir, xpath)

    h5createFile(ofile)
    writeSparseMatrix(x, ofile, name="sparse", column=TRUE, tenx=TRUE)
    nm <- .name_saver(x, ofile)

    list(
        `$schema`="hdf5_sparse_matrix/v1.json",
        path=xpath,
        is_child=child,
        `array` = .grab_array_type(x),
        hdf5_sparse_matrix = list(
            group = "sparse",
            format = "tenx_matrix",
            dimnames = nm
        )
    )
}

#' @importFrom DelayedArray is_sparse
.stage_any_matrix <- function(x, dir, path, child=FALSE) {
    if (is_sparse(x)) {
        .stage_sparse_matrix(x, dir, path, child=child)
    } else {
        .stage_array(x, dir, path, child=child)
    }
}

#' @export
#' @rdname stageArray
#' @importClassesFrom Matrix Matrix
setMethod("stageObject", "Matrix", function(x, dir, path, child=FALSE) .stage_any_matrix(x, dir, path, child=child))

#' @export
#' @rdname stageArray
setMethod("stageObject", "DelayedMatrix", function(x, dir, path, child=FALSE) .stage_delayed(x, dir, path, child = child, fallback = .stage_any_matrix))

#' @export
#' @rdname stageArray
recycleHdf5Files <- function(recycle) {
    prev <- staging.options$recycle.hdf5
    if (missing(recycle)) {
        prev
    } else {
        staging.options$recycle.hdf5 <- recycle
        invisible(prev)
    }
}

.link_or_copy <- function(from, to) {
    if (!file.link(from, to)) {
        if (!file.copy(from, to)) {
            stop("failed to link or copy '", from, "' to '", to, "'")
        }
    }
}

#' @importFrom BiocGenerics path
#' @importFrom rhdf5 h5read
.check_for_hdf5 <- function(x, dir, path, child) {
    if (!recycleHdf5Files()) {
        return(NULL)
    }

    if (is(x, "DelayedArray")) {
        x <- x@seed
    }
    if (is(x, "DelayedSetDimnames")) {
        x <- x@seed # Ignoring the names.
    }

    dense <- FALSE
    if (is(x, "HDF5ArraySeed")) {
        dense <- TRUE
    } else if (is(x, "H5SparseMatrixSeed")) {
        ;
    } else {
        return(NULL)
    }

    src <- path(x)

    if (!dense) {
        # Checking that it's in the 10X format, otherwise we bail.
        if (!is(x, "CSC_H5SparseMatrixSeed")) {
            return(NULL)
        }
        if (!is.null(x@subdata)) {
            return(NULL)
        }

        attempt <- try(h5read(src, paste0(x@group, "/shape")), silent=TRUE)
        if (is(attempt, "try-error") || !identical(as.vector(attempt), dim(x))) {
            return(NULL)
        }

        dir.create(file.path(dir, path), showWarnings=FALSE)
        dest <- paste0(path, "/matrix.h5")
        .link_or_copy(src, file.path(dir, dest))

        list(
            `$schema`="hdf5_sparse_matrix/v1.json",
            path=dest,
            is_child=child,
            `array` = .grab_array_type(x),
            hdf5_sparse_matrix = list(
                group = x@group,
                format = "tenx_matrix"
            )
        )

    } else {
        dir.create(file.path(dir, path), showWarnings=FALSE)
        dest <- paste0(path, "/array.h5")
        .link_or_copy(src, file.path(dir, dest))

        list(
            `$schema`="hdf5_dense_array/v1.json",
            path=dest,
            is_child=child,
            `array` = .grab_array_type(x),
            hdf5_dense_array= list(
                dataset=x@name
            )
        )
    }
}
