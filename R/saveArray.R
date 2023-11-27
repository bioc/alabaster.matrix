#' Save a multi-dimensional array to disk
#'
#' Save a high-dimensional array to its on-disk representations.
#'
#' @param x An integer, numeric, logical or character array.
#' Alternatively, any of the \linkS4class{denseMatrix} subclasses from the \pkg{Matrix} package.
#' @param path String containing the path to a directory in which to save \code{x}.
#' @param ... Further arguments, currently ignored.
#' 
#' @return
#' \code{x} is saved to \code{path} and \code{NULL} is invisibly returned.
#'
#' @seealso
#' \code{\link{readArray}}, to read the directory contents back into the R session.
#'
#' @author Aaron Lun
#' @examples
#' mat <- array(rpois(10000, 10), c(50, 20, 10))
#' dimnames(mat) <- list(
#'    paste0("GENE_", seq_len(nrow(mat))),
#'    letters[1:20],
#'    NULL
#' )
#'
#' dir <- tempfile()
#' saveObject(mat, dir)
#' list.files(dir)
#'
#' @name saveArray
#' @aliases 
#' stageObject,array-method
#' stageObject,Matrix-method
NULL

#' @import alabaster.base rhdf5
.save_array <- function(x, path, ...) {
    dir.create(path)
    fpath <- file.path(path, "array.h5")
    name <- "dense_array"

    # This needs to be wrapped up as writeHDF5Array needs to
    # take ownership of the file handle internally.
    local({
        fhandle <- H5Fcreate(fpath, "H5F_ACC_TRUNC")
        on.exit(H5Fclose(fhandle), add=TRUE, after=FALSE)

        ghandle <- H5Gcreate(fhandle, name)
        on.exit(H5Gclose(ghandle), add=TRUE, after=FALSE)

        h5_write_attribute(ghandle, "version", "1.0", scalar=TRUE)
        h5_write_attribute(ghandle, "type", array_type(x), scalar=TRUE)
        h5_write_attribute(ghandle, "transposed", 1L, scalar=TRUE)
    })

    transformed <- transformVectorForHdf5(x)
    writeHDF5Array(transformed$transformed, filepath=fpath, name="dense_array/data")

    # Reinitializing the handle for more downstream work.
    fhandle <- H5Fopen(fpath, "H5F_ACC_RDWR")
    on.exit(H5Fclose(fhandle), add=TRUE, after=FALSE)
    ghandle <- H5Gopen(fhandle, name)
    on.exit(H5Gclose(ghandle), add=TRUE, after=FALSE)

    if (!is.null(transformed$placeholder)) {
        dhandle <- H5Dopen(ghandle, "data")
        on.exit(H5Dclose(dhandle), add=TRUE, after=FALSE)
        h5_write_attribute(dhandle, missingPlaceholderName, transformed$placeholder, scalar=TRUE)
    }

    save_names(ghandle, x, transpose=TRUE)
    write(name, file=file.path(path, "OBJECT"))
    invisible(NULL)
}

#' @export
#' @rdname saveArray
setMethod("saveObject", "array", .save_array)

#' @export
#' @rdname saveArray
setMethod("saveObject", "denseMatrix", .save_array)

##############################
######### OLD STUFF ##########
##############################

#' @importFrom DelayedArray is_sparse
#' @importFrom rhdf5 h5createFile 
#' @importFrom HDF5Array writeHDF5Array
#' @importFrom alabaster.base transformVectorForHdf5 addMissingPlaceholderAttributeForHdf5
.stage_array <- function(x, dir, path, child=FALSE, .version=2) {
    dir.create(file.path(dir, path), showWarnings=FALSE)
    xpath <- paste0(path, "/array.h5")
    ofile <- file.path(dir, xpath)

    transformed <- transformVectorForHdf5(x)

    h5createFile(ofile)
    writeHDF5Array(transformed$transformed, filepath=ofile, name="data")
    nm <- .name_saver(x, ofile)

    if (!is.null(transformed$placeholder)) {
        addMissingPlaceholderAttributeForHdf5(ofile, "data", transformed$placeholder)
    }

    list(
        `$schema` = "hdf5_dense_array/v1.json",
        path = xpath,
        is_child = child,
        `array` = .grab_array_type(x),
        hdf5_dense_array = list(
            dataset = "data",
            dimnames = nm,
            version = 2
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
    xpath <- paste0(path, "/delayed.h5")
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
setMethod("stageObject", "DelayedArray", function(x, dir, path, child=FALSE) .stage_delayed(x, dir, path, child = child, fallback = .stage_array))

#' @importFrom rhdf5 h5createFile h5createGroup
#' @importFrom HDF5Array writeHDF5Array
.stage_sparse_matrix <- function(x, dir, path, child=FALSE) {
    dir.create(file.path(dir, path), showWarnings=FALSE)
    xpath <- paste0(path, "/matrix.h5")
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
            dimnames = nm,
            version = 2 
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
#' @importClassesFrom Matrix Matrix
setMethod("stageObject", "Matrix", function(x, dir, path, child=FALSE) .stage_any_matrix(x, dir, path, child=child))

#' @export
setMethod("stageObject", "DelayedMatrix", function(x, dir, path, child=FALSE) .stage_delayed(x, dir, path, child = child, fallback = .stage_any_matrix))

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
                format = "tenx_matrix",
                version = 2
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
                dataset=x@name,
                version = 2
            )
        )
    }
}
