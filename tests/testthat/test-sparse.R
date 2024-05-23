# This tests the sparse readers/writers.
# library(testthat); library(alabaster.matrix); source("test-sparse.R")

library(Matrix)
library(SparseArray)

test_that("writing to a sparse matrix works as expected", {
    for (i in 1:3) {
        x <- rsparsematrix(100, 20, 0.5)
        if (i == 2) {
            x <- DelayedArray(x) * 1 # force use of the block method.
        } else if (i == 3) {
            x <- BiocGenerics::cbind(DelayedArray(x), DelayedArray(rsparsematrix(100, 10, 0.5)))
        } else if (i == 4) {
            x <- as(x, "SVT_SparseMatrix")
        } else if (i == 5) {
            x <- BiocGenerics::cbind(DelayedArray(as(x, "SVT_SparseMatrix")), DelayedArray(as(rsparsematrix(100, 10, 0.5), "SVT_SparseMatrix")))
        }

        tmp <- tempfile(fileext=".h5")
        writeSparseMatrix(x, tmp, "csc_matrix")
        if (i <= 3) { # TODO: SVT matrices don't work with read_sparse_block(), for some reason.
            writeSparseMatrix(x, tmp, "csr_matrix", column=FALSE)
        }
        writeSparseMatrix(x, tmp, "tenx_matrix", tenx=TRUE)

        library(HDF5Array)
        expect_identical(unname(as.matrix(H5SparseMatrix(tmp, "csc_matrix"))), unname(as.matrix(x)))
        if (i <= 3) {
            expect_identical(unname(as.matrix(H5SparseMatrix(tmp, "csr_matrix"))), unname(as.matrix(x)))
        }
        expect_identical(unname(as.matrix(H5SparseMatrix(tmp, "tenx_matrix"))), unname(as.matrix(x)))
    }
})

test_that("writing to a sparse matrix works with tiny chunks", {
    for (i in 1:2) {
        x <- rsparsematrix(100, 20, 0.5)
        if (i == 2) {
            x <- DelayedArray(x) * 1 # force use of the block method.
        } else if (i == 3) {
            x <- as(x, "SVT_SparseMatrix")
        }

        tmp <- tempfile(fileext=".h5")
        setAutoBlockSize(max(dim(x))*8)
        writeSparseMatrix(x, tmp, "csc_matrix")
        if (i <= 2) { # TODO: uncomment when SVT matrices work with read_sparse_block().
            writeSparseMatrix(x, tmp, "csr_matrix", column=FALSE)
        }

        library(HDF5Array)
        expect_identical(unname(as.matrix(H5SparseMatrix(tmp, "csc_matrix"))), unname(as.matrix(x)))
        if (i <= 2) {
            expect_identical(unname(as.matrix(H5SparseMatrix(tmp, "csr_matrix"))), unname(as.matrix(x)))
        }
    }

    setAutoBlockSize()
})

get_type <- function(tmp, path) {
    fhandle <- H5Fopen(tmp)
    on.exit(H5Fclose(fhandle))
    dhandle <- H5Dopen(fhandle, path)
    on.exit(H5Dclose(dhandle), add=TRUE, after=FALSE)
    .Call("_getDatatypeName", H5Dget_type(dhandle), PACKAGE = "rhdf5")
}

test_that("writing to a sparse matrix works with guessed type", {
    set.seed(1000)
    for (i in 1:3) {
        core <- function() round(rsparsematrix(100, 20, 0.5))
        if (i == 1) {
            FUN <- function(f) f(core())
        } else if (i == 2) {
            FUN <- function(f) DelayedArray(f(core())) * 1 # force use of the block method.
        } else if (i == 3) {
            FUN <- function(f) BiocGenerics::cbind(DelayedArray(f(core())), DelayedArray(f(round(rsparsematrix(100, 10, 0.5)))))
        } else if (i == 4) {
            FUN <- function(f) as(f(core()), "SVT_SparseMatrix")
        }

        # Signed:
        tmp <- tempfile(fileext=".h5")
        y <- FUN(identity)
        writeSparseMatrix(y, tmp, "csc_matrix")
        expect_match(get_type(tmp, "csc_matrix/data"), "I32")
        expect_equal(unname(as.matrix(H5SparseMatrix(tmp, "csc_matrix"))), unname(as.matrix(y)))

        # Unsigned short
        tmp <- tempfile(fileext=".h5")
        y <- FUN(function(x) abs(x * 1000))
        writeSparseMatrix(y, tmp, "csc_matrix")
        expect_match(get_type(tmp, "csc_matrix/data"), "U16")
        expect_equal(unname(as.matrix(H5SparseMatrix(tmp, "csc_matrix"))), unname(as.matrix(y)))

        # Unsigned word
        tmp <- tempfile(fileext=".h5")
        y <- FUN(function(x) abs(x * 1000000))
        writeSparseMatrix(y, tmp, "csc_matrix")
        expect_match(get_type(tmp, "csc_matrix/data"), "I32")
        expect_equal(unname(as.matrix(H5SparseMatrix(tmp, "csc_matrix"))), unname(as.matrix(y)))

        # Unsigned long
        tmp <- tempfile(fileext=".h5")
        y <- FUN(function(x) abs(x * 10000000000))
        writeSparseMatrix(y, tmp, "csc_matrix")
        expect_match(get_type(tmp, "csc_matrix/data"), "F64")
        expect_identical(unname(as.matrix(H5SparseMatrix(tmp, "csc_matrix"))), unname(as.matrix(y)))
    }
})

test_that("writing to a sparse matrix works with guessed index type for block method", {
    set.seed(1000)
    for (i in 1:2) {
        for (big in c(TRUE, FALSE)) {
            nr <- if (big) 100000 else 100

            x <- round(rsparsematrix(nr, 20, 0.5))
            if (i == 2) {
                x <- DelayedArray(x) * 1 # force use of the block method.
            } else if (i == 3) {
                x <- as(x, "SVT_SparseMatrix")
            }

            tmp <- tempfile(fileext=".h5")
            writeSparseMatrix(x, tmp, "csc_matrix")
            expected.type <- if (big) "U32" else "U16"
            expect_match(get_type(tmp, "csc_matrix/indices"), expected.type)
        }
    }
})

test_that("writing to a sparse matrix works with NA values", {
    library(HDF5Array)
    x <- rsparsematrix(100, 20, 0.5)

    # No NA's added.
    {
        tmp <- tempfile(fileext=".h5")
        writeSparseMatrix(x, tmp, "csc_matrix")
        expect_null(rhdf5::h5readAttributes(tmp, "csc_matrix/data")[["missing-value-placeholder"]])
    }

    # Double-precision mode.
    x@x[1] <- NA
    {
        tmp <- tempfile(fileext=".h5")
        writeSparseMatrix(x, tmp, "csc_matrix")
        expect_identical(unname(as.matrix(H5SparseMatrix(tmp, "csc_matrix"))), unname(as.matrix(x)))
        expect_identical(rhdf5::h5readAttributes(tmp, "csc_matrix/data")[["missing-value-placeholder"]], NA_real_)
    }

    # Integer mode.
    x <- round(x)
    {
        tmp <- tempfile(fileext=".h5")
        writeSparseMatrix(x, tmp, "csc_matrix")
        expect_equal(unname(as.matrix(H5SparseMatrix(tmp, "csc_matrix"))), unname(as.matrix(x)))
        expect_identical(rhdf5::h5readAttributes(tmp, "csc_matrix/data")[["missing-value-placeholder"]], NA_integer_)
    }
})
