# This tests the sparse readers/writers.
# library(testthat); library(alabaster.matrix); source("setup.R"); source("test-SparseMatrix.R")

library(Matrix)
library(SparseArray)
library(DelayedArray)

test_that("reading a sparse matrix works with different output types", {
    x <- rsparsematrix(109, 297, 0.5)

    tmp <- tempfile()
    saveObject(x, tmp)
    roundtrip <- readObject(tmp)
    expect_identical(BiocGenerics::path(roundtrip), tmp)
    expect_s4_class(roundtrip, "ReloadedArray")
    expect_true(is_sparse(roundtrip))
    expect_identical(as(roundtrip, "dgCMatrix"), x)

    # Trying with a logical matrix.
    x <- rsparsematrix(109, 297, 0.5) > 0
    tmp <- tempfile()
    saveObject(x, tmp)
    roundtrip2 <- readObject(tmp)
    expect_identical(as.matrix(roundtrip2), as.matrix(x)) # use dense matrix as 'x' still contains zeros in its structural non-zeros.

    # Trying with an integer matrix.
    x <- rsparsematrix(109, 297, 0.5) * 10
    x <- as(x, "SVT_SparseMatrix")
    type(x) <- "integer"
    tmp <- tempfile()
    saveObject(x, tmp)
    roundtrip3 <- readObject(tmp)
    expect_identical(as(roundtrip3, "dgCMatrix"), as(x, "dgCMatrix"))
})

test_that("writing to a sparse matrix works as expected for numeric data", {
    for (i in 1:4) {
        for (miss in c(TRUE, FALSE)) {
            x <- rsparsematrix(100, 20, 0.5)
            if (miss) {
                x[1,1] <- NA
            }

            if (i == 2) {
                x <- as(x, "RsparseMatrix")
            } else if (i == 3) {
                x <- as(x, "SVT_SparseMatrix")
            } else if (i == 4) {
                x <- DelayedArray(x) * 1 # force use of the block method.
            }

            tmp <- tempfile(fileext=".h5")
            saveObject(x, tmp)
            roundtrip <- readObject(tmp)
            expect_identical_without_names(as.matrix(roundtrip), as.matrix(x))
        }
    }
})

test_that("writing to a sparse matrix works as expected for logical data", {
    for (i in 1:4) {
        for (miss in c(TRUE, FALSE)) {
            x <- rsparsematrix(100, 20, 0.5) > 0
            if (miss) {
                x[1,1] <- NA
            }

            if (i == 2) {
                x <- as(x, "RsparseMatrix")
            } else if (i == 3) {
                x <- as(x, "SVT_SparseMatrix")
            } else if (i == 4) {
                x <- DelayedArray(x) & TRUE # force use of the block method.
            }

            tmp <- tempfile(fileext=".h5")
            saveObject(x, tmp)
            roundtrip <- readObject(tmp)
            expect_identical_without_names(as.matrix(roundtrip), as.matrix(x))
        }
    }
})

test_that("writing to a sparse matrix works as expected for integer data", {
    for (i in 1:2) {
        for (miss in c(TRUE, FALSE)) {
            x <- round(rsparsematrix(100, 20, 0.5) * 10)
            if (miss) {
                x[1,1] <- NA
            }

            if (i == 1) {
                x <- as(x, "SVT_SparseMatrix")
            } else if (i == 2) {
                x <- DelayedArray(x) * 1L # force use of the block method.
            }
            type(x) <- "integer"

            tmp <- tempfile(fileext=".h5")
            saveObject(x, tmp)
            roundtrip <- readObject(tmp)
            expect_identical_without_names(as.matrix(roundtrip), as.matrix(x))
        }
    }
})

test_that("depositing a large sparseMatrix vector works correctly", {
    # We generate 1 million non-zero elements, which should exceed
    # the size of the chunks used inside the HDF5 file.
    x <- rsparsematrix(10000, 500, 0.2)
    expect_true(length(x@x) > h5_guess_vector_chunks(length(x@x)))

    tmp <- tempfile(fileext=".h5")
    saveObject(x, tmp)
    roundtrip <- readObject(tmp)
    expect_identical(as(roundtrip, 'dgCMatrix'), x)

    # Now injecting an NA, which should force it to use chunk-wise processing.
    x@x[1] <- NA

    tmp <- tempfile(fileext=".h5")
    saveObject(x, tmp)
    roundtrip <- readSparseMatrix(tmp)
#    expect_identical(as(roundtrip, 'dgCMatrix'), x) # TODO: bug in HDF5Array
})

test_that("depositing small chunks works correctly", {
    x <- rsparsematrix(1000, 500, 0.2)

    # For DelayedArrays:
    y <- DelayedArray(x) * 1 # force block processing.

    tmp <- tempfile()
    local({
        old <- getAutoBlockSize()
        setAutoBlockSize(2000)
        on.exit(setAutoBlockSize(old))
        saveObject(y, tmp)
    })

    roundtrip <- readObject(tmp)
    expect_identical(as(roundtrip, 'dgCMatrix'), x)

    # For SVT_SparseMatrix objects.
    y <- as(x, "SVT_SparseMatrix")

    tmp <- tempfile()
    local({
        old <- getAutoBlockSize()
        setAutoBlockSize(2000)
        on.exit(setAutoBlockSize(old))
        saveObject(y, tmp)
    })

    roundtrip <- readObject(tmp)
    expect_identical(as(roundtrip, 'dgCMatrix'), x)
})

test_that("fallback to large integer types for indices works correctly", {
    for (i in 1:4) {
        x <- rsparsematrix(100000, 20, 0.001)
        x[100000,20] <- 99 # making sure there's a value at the bottom-right so that we check the index correctly.

        if (i == 2) {
            x <- Matrix::t(x)
            x <- as(x, "RsparseMatrix")
        } else if (i == 3) {
            x <- as(x, "SVT_SparseMatrix")
        } else if (i == 4) {
            x <- DelayedArray(x) * 1 # force use of the block method.
        }

        tmp <- tempfile(fileext=".h5")
        saveObject(x, tmp)
        roundtrip <- readObject(tmp)
        expect_identical_without_names(as.matrix(roundtrip), as.matrix(x))
    }
})

test_that("writing to a sparse matrix works with different integer types", {
    for (i in 1:2) { # 1:3) { # bug in H5SparseMatrixSeed for INT types greater than int32.
        if (i == 1) {
            bits <- 8
        } else if (i == 2) {
            bits <- 16
        } else {
            bits <- 32
        }

        core <- matrix(0, 100, 200)
        chosen <- sample(20000, 1000)
        core[chosen] <- floor(runif(1000, 0, 2^bits))

        tmp <- tempfile()
        y <- as(core, "dgCMatrix")
        saveObject(y, tmp)
        expect_identical(as.matrix(readObject(tmp)), core)

        core[chosen] <- core[chosen] - 2^(bits-1)
        y <- as(core, "dgCMatrix")
        tmp <- tempfile()
        saveObject(y, tmp)
        expect_identical(as.matrix(readObject(tmp)), core)
    }
})

test_that("writing a DelayedMatrixworks with row-major storage", {
    xc <- rsparsematrix(100, 20, 0.5)
    core <- as.matrix(xc)

    # Regular:
    tmp <- tempfile()
    saveObject(xc, tmp)
    expect_identical(rhdf5::h5readAttributes(file.path(tmp, "matrix.h5"), "compressed_sparse_matrix")$layout, "CSC")
    expect_identical(unname(as.matrix(readObject(tmp))), core)

    # Saving as CSC based on chunkdims:
    setClass("TestColMatrix", contains="dgCMatrix")
    setMethod("chunkdim", "TestColMatrix", function(x) c(nrow(x), 1L))
    y <- as(xc, "TestColMatrix")
    y <- DelayedArray(y) * 1 # force the use of block processing.

    tmp <- tempfile()
    saveObject(y, tmp)
    expect_identical(rhdf5::h5readAttributes(file.path(tmp, "matrix.h5"), "compressed_sparse_matrix")$layout, "CSC")
    expect_identical(as.matrix(readObject(tmp)), core)

    # Saving as CSR:
    xr <- as(xc, "RsparseMatrix")

    setClass("TestRowMatrix", contains="dgRMatrix")
    setMethod("chunkdim", "TestRowMatrix", function(x) c(1L, ncol(x)))
    y <- as(xr, "TestRowMatrix")
    y <- DelayedArray(y) * 1 # force the use of block processing.

    tmp <- tempfile()
    saveObject(y, tmp)
    expect_identical(rhdf5::h5readAttributes(file.path(tmp, "matrix.h5"), "compressed_sparse_matrix")$layout, "CSR")
    expect_identical(as.matrix(readObject(tmp)), core)
})

test_that("writing an SVT matrix works with empty or near-empty storage", {
    for (it in 1:3) {
        if (it == 1L) {
            # Near-empty, with small integers.
            x <- rsparsematrix(100, 20, 0.001) * 10
        } else if (it == 2L) {
            # Near-empty, with larger integers. This checks that the type
            # optimization doesn't get screwed up by Infs after performing
            # range() on an empty vector.
            x <- rsparsematrix(100, 20, 0.001) * 1e6
        } else {
            # Actually empty.
            x <- rsparsematrix(100, 20, 0)
        }
        y <- as(x, "SVT_SparseMatrix")
        core <- as.matrix(x)

        tmp <- tempfile()
        saveObject(y, tmp)
        expect_identical(as.matrix(readObject(tmp)), core)

        # Now integers:
        y2 <- y
        type(y2) <- "integer"
        core2 <- core
        storage.mode(core2) <- "integer"

        tmp <- tempfile()
        saveObject(y2, tmp)
        expect_identical(as.matrix(readObject(tmp)), core2)

        # Now booleans:
        y2 <- y
        type(y2) <- "logical"
        core2 <- core
        storage.mode(core2) <- "logical"

        tmp <- tempfile()
        saveObject(y2, tmp)
        expect_identical(as.matrix(readObject(tmp)), core2)
    }
})

test_that("writing a sparse matrix works when the entire matrix consists of NAs", {
    x <- rsparsematrix(100, 20, 0)
    x@x[] <- NA # all NA's.
    y <- as(x, "SVT_SparseMatrix")
    core <- as.matrix(x)

    tmp <- tempfile()
    saveObject(y, tmp)
    expect_identical(as.matrix(readObject(tmp)), core)

    # Now integers:
    y2 <- y
    type(y2) <- "integer"
    core2 <- core
    storage.mode(core2) <- "integer"

    tmp <- tempfile()
    saveObject(y2, tmp)
    expect_identical(as.matrix(readObject(tmp)), core2)

    # Now booleans:
    y2 <- y
    type(y2) <- "logical"
    core2 <- core
    storage.mode(core2) <- "logical"

    tmp <- tempfile()
    saveObject(y2, tmp)
    expect_identical(as.matrix(readObject(tmp)), core2)
})

test_that("saving sparse matrices works with dimnames", {
    x <- rsparsematrix(100, 20, 0.5)
    rownames(x) <- paste0("GENE_", seq_len(nrow(x)))
    colnames(x) <- head(LETTERS, 20)

    tmp <- tempfile()
    saveObject(x, tmp)
    roundtrip <- readObject(tmp)
    expect_identical(as.matrix(roundtrip), as.matrix(x))
})

test_that("saveObject diverts correctly with pristine sparse DelayedArrays", {
    x <- DelayedArray(Matrix::rsparsematrix(100, 20, 0.2))
    expect_true(isPristine(x))

    tmp <- tempfile()
    saveObject(x, tmp)
    expect_identical(as(readObject(tmp), "dgCMatrix"), x@seed)
})

test_that("reading sparse arrays work with non-default NA placeholders", {
    x <- rsparsematrix(100, 20, 0.5)
    first <- x@x[1]

    dir <- tempfile()
    saveObject(x, dir)

    library(rhdf5)
    local({ 
        fhandle <- H5Fopen(file.path(dir, "matrix.h5"))
        on.exit(H5Fclose(fhandle), add=TRUE, after=FALSE)
        ghandle <- H5Gopen(fhandle, "compressed_sparse_matrix")
        on.exit(H5Gclose(ghandle), add=TRUE, after=FALSE)
        dhandle <- H5Dopen(ghandle, "data")
        on.exit(H5Dclose(dhandle), add=TRUE, after=FALSE)
        alabaster.base::h5_write_attribute(dhandle, "missing-value-placeholder", first, scalar=TRUE)
    })

    arr2 <- readObject(dir)
    ref <- as.matrix(x)
    ref[ref==first] <- NA
    expect_identical(ref, as.matrix(arr2))
})
