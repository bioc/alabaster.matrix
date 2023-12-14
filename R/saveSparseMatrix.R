#' Save a sparse matrix to disk
#'
#' Save a sparse matrix to its on-disk representations.
#'
#' @param x A sparse matrix of some kind, typically from either the \pkg{Matrix} or \pkg{SparseArray} packages.
#' @param path String containing the path to a directory in which to save \code{x}.
#' @param ... Further arguments, currently ignored.
#' 
#' @return
#' \code{x} is saved to \code{path} and \code{NULL} is invisibly returned.
#'
#' @seealso
#' \code{\link{readSparseMatrix}}, to read the directory contents back into the R session.
#'
#' @author Aaron Lun
#' @examples
#' mat <- Matrix::rsparsematrix(100, 200, density=0.2)
#' rownames(mat) <- paste0("GENE_", seq_len(nrow(mat)))
#'
#' dir <- tempfile()
#' saveObject(mat, dir)
#' list.files(dir)
#'
#' @name saveSparseMatrix
NULL

.save_compressed_sparse_matrix <- function(x, path, ...) {
    dir.create(path)
    fpath <- file.path(path, "matrix.h5")
    name <- "compressed_sparse_matrix"

    fhandle <- H5Fcreate(fpath, "H5F_ACC_TRUNC")
    on.exit(H5Fclose(fhandle), add=TRUE, after=FALSE)

    ghandle <- H5Gcreate(fhandle, name)
    on.exit(H5Gclose(ghandle), add=TRUE, after=FALSE)

    h5_write_attribute(ghandle, "type", to_array_type(x), scalar=TRUE)
    h5_write_vector(ghandle, "shape", dim(x), type="H5T_NATIVE_UINT32")

    optimized <- optimize_storage(x)
    h5_write_sparse_matrix(x, handle=ghandle, details=optimized)

    if (!is.null(optimized$placeholder)) {
        dhandle <- H5Dopen(ghandle, "data")
        on.exit(H5Dclose(dhandle), add=TRUE, after=FALSE)
        h5_write_attribute(dhandle, missingPlaceholderName, optimized$placeholder, type=optimized$type, scalar=TRUE)
    }

    save_names(ghandle, x)
    saveObjectFile(path, name, list(compressed_sparse_matrix=list(version="1.0")))
    invisible(NULL)
}

#' @export
#' @rdname saveSparseMatrix
setMethod("saveObject", "sparseMatrix", function(x, path, ...) .save_compressed_sparse_matrix(x, path, ...))

#' @export
#' @rdname saveSparseMatrix
setMethod("saveObject", "SVT_SparseMatrix", function(x, path, ...) .save_compressed_sparse_matrix(x, path, ...))

##############################
######### INTERNALS ##########
##############################

setGeneric("h5_write_sparse_matrix", function(x, ...) standardGeneric("h5_write_sparse_matrix"))

.deposit_sparse_data_with_placeholder <- function(handle, nzd, type, placeholder) {
    N <- length(nzd)
    chunks <- h5_guess_vector_chunks(N)
    dhandle <- h5_create_vector(handle, "data", N, chunks=chunks, type=type)
    on.exit(H5Dclose(dhandle), add=TRUE, after=FALSE)

    shandle <- H5Screate_simple(N)
    on.exit(H5Sclose(shandle), add=TRUE, after=FALSE)
    mem_shandle <- H5Screate_simple(chunks)
    on.exit(H5Sclose(mem_shandle), add=TRUE, after=FALSE)

    pos <- 1L
    while (pos <= N) {
        size <- chunks
        end <- pos + size - 1L
        if (end > N) {
            size <- N - pos + 1L
            end <- N
            H5Sset_extent_simple(mem_shandle, size)
        }
        H5Sselect_hyperslab(shandle, "H5S_SELECT_SET", start=pos, count=size)

        current <- nzd[pos:end]
        current[is.missing(current)] <- placeholder
        H5Dwrite(dhandle, current, h5spaceMem=mem_shandle, h5spaceFile=shandle)
        pos <- end + 1L
    }
}

.choose_itype <- function(ilimit) {
    if (ilimit < 2^16) {
        "H5T_NATIVE_UINT16"
    } else {
        "H5T_NATIVE_UINT32"
    }
}

setMethod("h5_write_sparse_matrix", "CsparseMatrix", function(x, handle, details, ...) {
    if (!is.null(details$placeholder)) {
        .deposit_sparse_data_with_placeholder(handle, x@x, type=details$type, placeholder=details$placeholder)
    } else {
        h5_write_vector(handle, "data", x@x, type=details$type)
    }

    h5_write_vector(handle, "indices", x@i, type=.choose_itype(nrow(x)))
    h5_write_vector(handle, "indptr", x@p, type="H5T_NATIVE_UINT64")
    h5_write_attribute(handle, "layout", "CSC", scalar=TRUE)
})

setMethod("h5_write_sparse_matrix", "RsparseMatrix", function(x, handle, details, ...) {
    if (!is.null(details$placeholder)) {
        .deposit_sparse_data_with_placeholder(handle, x@x, type=details$type, placeholder=details$placeholder)
    } else {
        h5_write_vector(handle, "data", x@x, type=details$type)
    }

    h5_write_vector(handle, "indices", x@j, type=.choose_itype(ncol(x)))
    h5_write_vector(handle, "indptr", x@p, type="H5T_NATIVE_UINT64")
    h5_write_attribute(handle, "layout", "CSR", scalar=TRUE)
})

setMethod("h5_write_sparse_matrix", "SVT_SparseMatrix", function(x, handle, details, ...) {
    N <- details$size
    chunks <- h5_guess_vector_chunks(N)
    dhandle <- h5_create_vector(handle, "data", N, chunks=chunks, type=details$type)
    on.exit(H5Dclose(dhandle), add=TRUE, after=FALSE)

    ihandle <- h5_create_vector(handle, "indices", N, chunks=chunks, type=.choose_itype(nrow(x)))
    on.exit(H5Dclose(ihandle), add=TRUE, after=FALSE)

    shandle <- H5Screate_simple(N)
    on.exit(H5Sclose(shandle), add=TRUE, after=FALSE)
    mhandle <- H5Screate_simple(1)
    on.exit(H5Sclose(mhandle), add=TRUE, after=FALSE)

    SVT <- x@SVT
    if (is.null(SVT)) {
        column.counts <- integer(ncol(x))

    } else {
        start <- 1L
        cached <- 0
        last.cleared <- 0L
        chunksize <- getAutoBlockLength(type(x))

        column.counts <- integer(length(SVT))
        for (i in seq_along(SVT)) {
            y <- SVT[[i]]
            if (!is.null(y)) {
                column.counts[i] <- length(y[[1]])
            }
        }

        for (i in seq_len(ncol(x))) {
            cached <- cached + column.counts[i]

            if (cached >= chunksize || i == ncol(x)) {
                targets <- SVT[(last.cleared + 1L):i]
                all.i <- all.d <- vector("list", length(targets))
                for (j in seq_along(targets)) {
                    y <- targets[[j]]
                    if (is.null(y)) {
                        all.i[[j]] <- integer(0)
                        all.d[[j]] <- as(NULL, type(x))
                    } else {
                        all.i[[j]] <- y[[1]]
                        all.d[[j]] <- y[[2]]
                    }
                }

                all.i <- unlist(all.i)
                all.d <- unlist(all.d)
                stopifnot(cached == length(all.d))

                if (!is.null(details$placeholder)) {
                    all.d[is.missing(all.d)] <- details$placeholder
                }

                H5Sselect_hyperslab(shandle, "H5S_SELECT_SET", start=start, count=cached)
                H5Sset_extent_simple(mhandle, cached)
                H5Dwrite(dhandle, all.d, h5spaceMem=mhandle, h5spaceFile=shandle)
                H5Dwrite(ihandle, all.i, h5spaceMem=mhandle, h5spaceFile=shandle)

                last.cleared <- i
                start <- start + cached
                cached <- 0
            }
        }
    }

    h5_write_vector(handle, "indptr", c(0, cumsum(column.counts)), type="H5T_NATIVE_UINT64")
    h5_write_attribute(handle, "layout", "CSC", scalar=TRUE)
})

setMethod("h5_write_sparse_matrix", "ANY", function(x, handle, details, ...) {
    cd <- chunkdim(x)
    if (is.null(cd) || !isFALSE(ncol(x) / cd[2] >= nrow(x) / cd[1])) { # i.e., number of chunks per row >= number of chunks per column, indicating it's easier to extract by column.
        layout <- "CSC"
        ilimit <- nrow(x)
        grid <- colAutoGrid(x)
    } else {
        layout <- "CSR"
        ilimit <- ncol(x)
        grid <- rowAutoGrid(x)
    }

    out <- vector("list", length(grid))

    N <- details$size
    chunks <- h5_guess_vector_chunks(N)
    dhandle <- h5_create_vector(handle, "data", N, type=details$type, chunks=chunks)
    on.exit(H5Dclose(dhandle), add=TRUE, after=FALSE)

    ihandle <- h5_create_vector(handle, "indices", N, type=.choose_itype(ilimit), chunks=chunks)
    on.exit(H5Dclose(ihandle), add=TRUE, after=FALSE)

    shandle <- H5Screate_simple(N)
    on.exit(H5Sclose(shandle), add=TRUE, after=FALSE)
    mhandle <- H5Screate_simple(1)
    on.exit(H5Sclose(mhandle), add=TRUE, after=FALSE)

    start <- 1L
    pointers <- list(0)
    for (i in seq_along(grid)) {
        block <- read_sparse_block(x, grid[[i]])

        nzdex <- nzindex(block)
        if (layout == "CSC") {
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

        if (!is.null(details$placeholder)) {
            v[is.missing(v)] <- details$placeholder
        }

        H5Sselect_hyperslab(shandle, "H5S_SELECT_SET", start=start, count=length(v))
        H5Sset_extent_simple(mhandle, length(v))
        H5Dwrite(dhandle, v, h5spaceMem=mhandle, h5spaceFile=shandle)
        H5Dwrite(ihandle, secondary - 1L, h5spaceMem=mhandle, h5spaceFile=shandle)
        pointers <- c(pointers, list(tabulate(primary, ndim)))
        start <- start + length(v)
    }

    h5_write_vector(handle, "indptr", cumsum(unlist(pointers)), type="H5T_NATIVE_UINT64")
    h5_write_attribute(handle, "layout", layout, scalar=TRUE)
})
