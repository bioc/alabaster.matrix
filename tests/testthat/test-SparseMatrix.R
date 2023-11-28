# This tests the sparse readers/writers.
# library(testthat); library(alabaster.matrix); source("test-SparseMatrix.R")

library(Matrix)
library(SparseArray)

test_that("writing to a sparse matrix works as expected for numeric data", {
    for (i in 1:3) {
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
            expect_identical(as.matrix(roundtrip), as.matrix(x))
        }
    }
})

test_that("writing to a sparse matrix works as expected for logical data", {
    for (i in 1:3) {
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
            expect_identical(as.matrix(roundtrip), as.matrix(x))
        }
    }
})

test_that("writing to a sparse matrix works as expected for integer data", {
    for (i in 1:1) {
        for (miss in c(TRUE, FALSE)) {
            x <- round(rsparsematrix(100, 20, 0.5) * 10)
            if (miss) {
                x[1,1] <- NA
            }

            if (i == 1) {
                x <- as(x, "SVT_SparseMatrix")
                type(x) <- "integer"
            } else if (i == 2) {
                x <- DelayedArray(x) * 1L # force use of the block method.
            }

            tmp <- tempfile(fileext=".h5")
            saveObject(x, tmp)
            roundtrip <- readObject(tmp)
            expect_identical(as.matrix(roundtrip), as.matrix(x))
        }
    }
})

test_that("depositing a large sparseMatrix vector works correctly", {
    # We generate 1 million non-zero elements, which should exceed
    # the size of the chunks used inside the HDF5 file.
    x <- rsparsematrix(10000, 500, 0.2)
    expect_true(length(x@x) > h5_guess_chunk_size(length(x@x)))

    tmp <- tempfile(fileext=".h5")
    saveObject(x, tmp)
    roundtrip <- readObject(tmp)
    expect_identical(as(roundtrip, 'dgCMatrix'), x)

    # Now injecting an NA.
    x@x[1] <- NA

    tmp <- tempfile(fileext=".h5")
    saveObject(x, tmp)
    roundtrip <- readSparseMatrix(tmp)
    expect_identical(as(roundtrip, 'dgCMatrix'), x)
})

test_that("fallback to large integer types for indices works correctly", {
    for (i in 1:3) {
        x <- rsparsematrix(100000, 20, 0.001)
        x[100000,20] <- 99 # making sure there's a value at the bottom-right so that we check the index correctly.

        if (i == 2) {
            x <- t(x)
            x <- as(x, "RsparseMatrix")
        } else if (i == 3) {
            x <- as(x, "SVT_SparseMatrix")
        } else if (i == 4) {
            x <- DelayedArray(x) * 1 # force use of the block method.
        }

        tmp <- tempfile(fileext=".h5")
        saveObject(x, tmp)
        roundtrip <- readObject(tmp)
        expect_identical(as.matrix(roundtrip), as.matrix(x))
    }
})

test_that("writing to a sparse matrix works with different integer types", {
    for (i in 1:3) {
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
