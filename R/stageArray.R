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
.stage_array <- function(x, dir, path, child=FALSE) {
    dir.create(file.path(dir, path), showWarnings=FALSE)
    xpath <- file.path(path, "array.h5")
    ofile <- file.path(dir, xpath)

    h5createFile(ofile)
    writeHDF5Array(x, filepath=ofile, name="data")
    nm <- .name_saver(x, ofile)

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
                h5write(current, file=ofile, name=file.path("names", i))
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
#' @importFrom DelayedArray DelayedArray seed isPristine
.stage_delayed <- function(x, dir, path, child, fallback) {
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

preserve.env <- new.env()
preserve.env$preserve <- FALSE

#' @export
#' @rdname stageArray
preserveDelayedOperations <- function(preserve) {
    prev <- preserve.env$preserve
    if (missing(preserve)) {
        prev
    } else {
        preserve.env$preserve <- preserve
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
        `$schema`="hdf5sparse_matrix/v1.json",
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
#' @rdname stageMatrix
#' @importClassesFrom Matrix Matrix
setMethod("stageObject", "Matrix", function(x, dir, path, child=FALSE) .stage_matrix(x, dir, path, child=child))

#' @export
#' @rdname stageMatrix
setMethod("stageObject", "DelayedMatrix", function(x, dir, path, child=FALSE) .stage_delayed(x, dir, path, child = child, fallback = .stage_any_matrix))
