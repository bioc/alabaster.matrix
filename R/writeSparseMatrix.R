#' Write a sparse matrix
#'
#' Writes a sparse matrix to file in a compressed sparse format.
#'
#' @param x A sparse matrix of some sort.
#' This includes sparse \linkS4class{DelayedMatrix} objects.
#' @param file String containing a path to the HDF5 file.
#' The file is created if it is not already present.
#' @param name String containing the name of the group to store \code{x}.
#' @param chunk Integer scalar specifying the chunk size for the indices and values.
#' @param column Logical scalar indicating whether to store as compressed sparse column format.
#' @param tenx Logical scalar indicating whether to use the 10X compressed sparse column format.
#' @param guess.integer Logical scalar specifying whether to guess an appropriate integer type from \code{x}.
#'
#' @details
#' This writes a sparse matrix to file in various formats:
#' \itemize{
#' \item \code{column=TRUE} and \code{tenx=FALSE} uses H5AD's \code{csr_matrix} format.
#' \item \code{column=FALSE} and \code{tenx=FALSE} uses H5AD's \code{csc_matrix} format.
#' \item \code{tenx=TRUE} uses 10X Genomics' HDF5 matrix format.
#' }
#' For the first two formats, the apparent transposition is deliberate, because columns in R are interpreted as rows in H5AD.
#' This allows us to retain consistency the interpretation of samples (columns in R, rows in H5AD) and features (vice versa).
#' Constructors for classes like \linkS4class{H5SparseMatrix} will automatically transpose so no extra work is required.
#'
#' If \code{guess.integer=TRUE}, we attempt to save \code{x}'s values into the smallest type that will accommodate all of its values. 
#' If \code{x} only contains unsigned integers, we will attempt to save either 8-, 16- or 32-bit unsigned integers.
#' If \code{x} contains signed integers, we will fall back to 32-bit signed integers.
#' For all other values, we will fall back to double-precision floating point values.
#'
#' We attempt to save \code{x}'s indices to unsigned 16-bit integers if the relevant dimension of \code{x} is small enough.
#' Otherwise we will save it as an unsigned 32-bit integer.
#'
#' @return
#' A \code{NULL} invisibly.
#' The contents of \code{x} are written to \code{name} in \code{file}.
#' 
#' @author Aaron Lun
#'
#' @examples
#' library(Matrix)
#' x <- rsparsematrix(100, 20, 0.5)
#' tmp <- tempfile(fileext=".h5")
#' writeSparseMatrix(x, tmp, "csc_matrix")
#' writeSparseMatrix(x, tmp, "csr_matrix", column=FALSE)
#' writeSparseMatrix(x, tmp, "tenx_matrix", tenx = TRUE)
#' 
#' rhdf5::h5ls(tmp)
#' library(HDF5Array)
#' H5SparseMatrix(tmp, "csc_matrix")
#' H5SparseMatrix(tmp, "csr_matrix")
#' H5SparseMatrix(tmp, "tenx_matrix")
#' @export
#' @importFrom rhdf5 h5createGroup h5createFile
writeSparseMatrix <- function(x, file, name, chunk=10000, column=TRUE, tenx=FALSE, guess.integer=TRUE) {
    if (!file.exists(file)) {
        h5createFile(file)
    }
    .write_CS_matrix(file, name, x, chunk_dim=chunk, by_column=column, use_tenx=tenx, guess_type=guess.integer)
    invisible(NULL)
}

#' @import methods
#' @importFrom rhdf5 h5write h5createGroup H5Gopen H5Gopen H5Fopen H5Fclose 
#' h5writeAttribute h5createDataset h5writeDataset H5Sunlimited H5Gclose
#' @importFrom DelayedArray colAutoGrid rowAutoGrid read_sparse_block type
#' @importClassesFrom Matrix dgCMatrix
.write_CS_matrix <- function(file, name, mat, chunk_dim = 10000, by_column=TRUE, use_tenx=FALSE, guess_type=TRUE) {
    handle <- H5Fopen(file)
    on.exit(H5Fclose(handle))
    h5createGroup(handle, name)

    if (use_tenx) {
        by_column <- TRUE # must be column major.
        h5writeDataset(dim(mat), handle, paste0(name, "/shape"))
    } else {
        if (by_column) {
            format <- "csr_matrix"
        } else {
            format <- "csc_matrix"
        }

        ghandle <- H5Gopen(handle, name)
        on.exit(H5Gclose(ghandle), add = TRUE, after = FALSE)

        h5writeAttribute(format, ghandle, "encoding-type", variableLengthString=TRUE, encoding="UTF8", asScalar=TRUE)
        h5writeAttribute("0.1.0", ghandle, "encoding-version", variableLengthString=TRUE, encoding="UTF8", asScalar=TRUE)
        h5writeAttribute(rev(dim(mat)), ghandle, "shape")
    }

    # Scan through and guess the best types.
    details <- .extract_sparse_details(mat)
    any.neg <- details$negative
    extreme <- details$extreme
    any.nonint <- details$non.integer
    count <- details$count
    has.missing <- details$missing

    if (any.nonint) {
        type <- "DOUBLE"
    } else if (extreme >= 2^31) {
        # Don't attempt to use unsigned 32-bit ints; values between 2^31
        # and 2^32 will fail to be converted to R's signed integers.
        type <- "DOUBLE"
    } else if (has.missing) {
        # Don't bother choosing a smaller type, we need to encode NA_integer_ as -2^31.
        type <- "INT32" 
    } else if (any.neg) {
        type <- "INT32"
    } else {
        # We don't try to save it as a byte for now, because that gets
        # interpreted by rhdf5 as a raw vector... not helpful.
        if (extreme < 2^16) {
            type <- "UINT16"
        } else {
            # Again, don't attempt to use unsigned 32-bit ints.
            type <- "INT32"
        }
    }

    htype <- paste0("H5T_NATIVE_", type)

    if (by_column) {
        max.index <- nrow(mat)
    } else {
        max.index <- ncol(mat)
    }

    if (max.index < 2^16) {
        itype <- "H5T_NATIVE_UINT16"
    } else {
        itype <- "H5T_NATIVE_UINT32"
    }

    chunk_dim <- min(details$count, chunk_dim)

    h5createDataset(
        handle,
        paste0(name, "/data"),
        dims = details$count,
        H5type = htype,
        chunk = chunk_dim
    )

    is_integer <- type != "DOUBLE"
    transformer <- identity
    if (has.missing) {
        if (is_integer) {
            transformer <- as.integer
        } else {
            transformer <- as.double
        }
        addMissingPlaceholderAttributeForHdf5(file, paste0(name, "/data"), transformer(NA))
    }

    h5createDataset(
        handle,
        paste0(name, "/indices"),
        dims = details$count,
        H5type = itype,
        chunk = chunk_dim
    )

    if (by_column) {
        p <- .dump_column_sparse_matrix(
            mat, 
            handle, 
            data.path=paste0(name, "/data"), 
            index.path=paste0(name, "/indices"), 
            start=NULL,
            transformer=transformer
        )
    } else {
        p <- .dump_row_sparse_matrix(
            mat, 
            handle, 
            data.path=paste0(name, "/data"), 
            index.path=paste0(name, "/indices"), 
            start=NULL,
            transformer=transformer
        )
    }
    indptrs <- c(0, cumsum(as.double(p)))

    iname <- paste0(name, "/indptr")
    h5createDataset(
        handle,
        iname,
        dims = length(indptrs),
        H5type = "H5T_NATIVE_UINT64"
    )

    h5writeDataset(indptrs, handle, iname)
}

####################################################

setGeneric(".extract_sparse_details", function(x) standardGeneric(".extract_sparse_details"))

.check_for_missing_value <- function(x) {
    if (anyNA(x)) {
        if (is.double(x)) {
            if (sum(is.na(x)) > sum(is.nan(x))) {
                return(TRUE)
            }
        } else {
            return(TRUE)
        }
    }
    return(FALSE)
}

#' @importClassesFrom Matrix dsparseMatrix
#' @importFrom DelayedArray getAutoBlockLength type
setMethod(".extract_sparse_details", "dsparseMatrix", function(x) {
    chunksize <- getAutoBlockLength(type(x))

    any.nonint <- FALSE
    i <- 1L
    while (i <= length(x@x)) {
        end <- i + chunksize - 1L
        current <- x@x[i:min(end, length(x@x))]
        if (any(current%%1!=0, na.rm=TRUE)) {
            any.nonint <- TRUE
            break
        } 
        i <- end + 1L
    }

    has.missing <- FALSE
    i <- 1L
    while (i <= length(x@x)) {
        end <- i + chunksize - 1L
        current <- x@x[i:min(end, length(x@x))]
        if (.check_for_missing_value(current)) {
            has.missing <- TRUE
            break
        }
        i <- end + 1L
    }

    limits <- range(x@x, na.rm=TRUE)
    list(
        missing = has.missing,
        negative = limits[1] < 0, 
        extreme = max(abs(limits), na.rm=TRUE),
        non.integer = any.nonint,
        count = length(x@x)
    )
})

#' @importClassesFrom SparseArray SVT_SparseMatrix
setMethod(".extract_sparse_details", "SVT_SparseMatrix", function(x) {
    limits <- range(unlist(lapply(x@SVT, function(x) range(x[[2]], na.rm=TRUE))), na.rm=TRUE)
    non.int <- any(vapply(x@SVT, function(x) any(x[[2]]%%1 != 0, na.rm=TRUE), TRUE), na.rm=TRUE)
    has.missing <- any(vapply(x@SVT, function(x) .check_for_missing_value(x[[2]]), TRUE), na.rm=TRUE)
    count <- sum(vapply(x@SVT, function(x) length(x[[1]]), 0L), na.rm=TRUE)

    list(
        missing = has.missing,
        negative = limits[1] < 0, 
        extreme = max(abs(limits), na.rm=TRUE),
        non.integer = non.int,
        count = count
    )
})

#' @importClassesFrom DelayedArray DelayedSetDimnames
setMethod(".extract_sparse_details", "DelayedSetDimnames", function(x) .extract_sparse_details(x@seed))

.combine_extracted_details <- function(collected) {
    list(
        missing = any(vapply(collected, function(y) y$missing, TRUE), na.rm=TRUE),
        negative = any(vapply(collected, function(y) y$negative, TRUE), na.rm=TRUE),
        extreme = max(unlist(lapply(collected, function(y) y$extreme)), na.rm=TRUE),
        non.integer = any(vapply(collected, function(y) y$non.integer, TRUE), na.rm=TRUE),
        count = sum(unlist(lapply(collected, function(y) y$count)), na.rm=TRUE)
    )
}

#' @importClassesFrom DelayedArray DelayedAbind
setMethod(".extract_sparse_details", "DelayedAbind", function(x) {
    collected <- lapply(x@seeds, .extract_sparse_details)
    .combine_extracted_details(collected)
})

#' @importClassesFrom DelayedArray DelayedMatrix
setMethod(".extract_sparse_details", "DelayedMatrix", function(x) .extract_sparse_details(x@seed))

#' @importFrom DelayedArray nzdata
.extract_sparse_details_fragment <- function(sparse) {
    vals <- nzdata(sparse)
    any.negative <- any(vals < 0, na.rm=TRUE)
    any.nonint <- any(vals != round(vals), na.rm=TRUE)
    extreme.val <- max(abs(vals), na.rm=TRUE)
    has.missing <- .check_for_missing_value(vals)

    list(
        missing = has.missing,
        negative = any.negative, 
        non.integer = any.nonint, 
        extreme = extreme.val, 
        count = length(vals)
    )
}

#' @importFrom DelayedArray blockApply
setMethod(".extract_sparse_details", "ANY", function(x) {
    collected <- blockApply(x, FUN = .extract_sparse_details_fragment, as.sparse = TRUE) # any grid works here.
    .combine_extracted_details(collected)
})

####################################################

setGeneric(".dump_column_sparse_matrix", function(x, handle, index.path, data.path, start, transformer) {
    standardGeneric(".dump_column_sparse_matrix")
})

#' @importFrom rhdf5 h5writeDataset
#' @importClassesFrom Matrix dgCMatrix
setMethod(".dump_column_sparse_matrix", "dgCMatrix", function(x, handle, index.path, data.path, start, transformer) {
    if (is.null(start)) {
        h5writeDataset(transformer(x@x), handle, data.path)
        h5writeDataset(x@i, handle, index.path)
    } else {
        index <- list(as.double(start) + seq_along(x@x))
        h5writeDataset(transformer(x@x), handle, data.path, index=index)
        h5writeDataset(x@i, handle, index.path, index=index)
    }
    diff(x@p)
})

.sanitize_start <- function(start) {
    if (is.null(start)) {
        0
    } else {
        as.double(start) # avoid integer overflow
    }
}

#' @importFrom rhdf5 h5writeDataset
#' @importClassesFrom SparseArray SVT_SparseMatrix
#' @importFrom DelayedArray getAutoBlockSize type 
setMethod(".dump_column_sparse_matrix", "SVT_SparseMatrix", function(x, handle, index.path, data.path, start, transformer) {
    # Processing things in chunks to reduce the number of HDF5 calls.
    chunksize <- min(getAutoBlockLength("integer"), getAutoBlockLength(type(x)))
    column.counts <- vapply(x@SVT, function(y) length(y[[1]]), 0L)

    start <- .sanitize_start(start)
    cached <- 0
    last.cleared <- 0L

    for (i in seq_len(ncol(x))) {
        cached <- cached + column.counts[i]

        if (cached >= chunksize || i == ncol(x)) {
            targets <- x@SVT[(last.cleared + 1L):i]
            all.i <- unlist(lapply(targets, function(y) y[[1]]))
            all.d <- unlist(lapply(targets, function(y) y[[2]]))

            index <- list(start + seq_along(all.i))
            h5writeDataset(all.i, handle, index.path, index = index)
            h5writeDataset(transformer(all.d), handle, data.path, index = index)

            last.cleared <- i
            cached <- 0
            start <- start + length(all.i)
        }
    }

    column.counts
})

#' @importClassesFrom DelayedArray DelayedSetDimnames
setMethod(".dump_column_sparse_matrix", "DelayedSetDimnames", function(x, handle, index.path, data.path, start, transformer) {
    .dump_column_sparse_matrix(x@seed, handle, index.path, data.path, start, transformer)
})

#' @importClassesFrom DelayedArray DelayedMatrix
setMethod(".dump_column_sparse_matrix", "DelayedMatrix", function(x, handle, index.path, data.path, start, transformer) {
    .dump_column_sparse_matrix(x@seed, handle, index.path, data.path, start, transformer)
})

#' @importClassesFrom DelayedArray DelayedAbind
setMethod(".dump_column_sparse_matrix", "DelayedAbind", function(x, handle, index.path, data.path, start, transformer) {
    if (x@along != 2L) {
        return(callNextMethod()) # goes to the ANY method.
    }

    start <- .sanitize_start(start)
    collected <- vector("list", length(x@seeds))
    for (s in seq_along(x@seeds)) {
        current <- .dump_column_sparse_matrix(x@seeds[[s]], handle, index.path, data.path, start=start, transformer=transformer)
        collected[[s]] <- current
        start <- start + sum(as.double(current))
    }

    unlist(collected)
})

#' @importFrom DelayedArray colAutoGrid read_sparse_block
setMethod(".dump_column_sparse_matrix", "ANY", function(x, handle, index.path, data.path, start, transformer) {
    start <- .sanitize_start(start)
    grid <- colAutoGrid(x)
    out <- vector("list", length(grid))

    for (i in seq_along(grid)) {
        block <- read_sparse_block(x, grid[[i]])
        cout <- .blockwise_sparse_writer(
            block, 
            start, 
            file=handle, 
            transformer=transformer,
            index.path=index.path, 
            data.path=data.path, 
            by_column=TRUE
        )
        
        out[[i]] <- cout$number
        start <- cout$last
    }

    unlist(out)
})

#' @importFrom DelayedArray rowAutoGrid read_sparse_block
.dump_row_sparse_matrix <- function(x, handle, index.path, data.path, start, transformer) {
    start <- .sanitize_start(start)
    grid <- rowAutoGrid(x)
    out <- vector("list", length(grid))

    for (i in seq_along(grid)) {
        block <- read_sparse_block(x, grid[[i]])
        cout <- .blockwise_sparse_writer(
            block, 
            start, 
            transformer=transformer, 
            file=handle, 
            index.path=index.path, 
            data.path=data.path, 
            by_column=FALSE
        )

        out[[i]] <- cout$number
        start <- cout$last
    }

    unlist(out)
}

#' @importFrom DelayedArray nzindex nzdata
#' @importFrom rhdf5 h5writeDataset
.blockwise_sparse_writer <- function(block, last, transformer, file, index.path, data.path, by_column) {
    nzdex <- nzindex(block)
    if (by_column) {
        primary <- nzdex[, 2]
        secondary <- nzdex[, 1]
        ndim <- ncol(block)
    } else {
        primary <- nzdex[, 1]
        secondary <- nzdex[, 2]
        ndim <- nrow(block)
    }
    v <- nzdata(block)

    o <- order(primary, secondary)
    primary <- primary[o]
    secondary <- secondary[o]
    v <- v[o]

    newlast <- last + length(primary)
    index <- list(last + seq_along(primary))
    h5writeDataset(secondary - 1L, file, index.path, index = index)
    h5writeDataset(transformer(v), file, data.path, index = index)

    list(number=tabulate(primary, ndim), last=newlast)
}
