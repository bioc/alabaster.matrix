# This checks that the storage optimization works as expected.
# library(testthat); library(alabaster.matrix); source("test-optimize_storage.R", encoding="UTF-8")

library(DelayedArray)
library(rhdf5)

test_that("storage optimization works for integers, no missing values", {
    old <- getAutoBlockSize()
    setAutoBlockSize(400)
    on.exit(setAutoBlockSize(old))

    for (i in 1:2) {
        fun <- identity
        if (i == 1L) {
            fun <- DelayedArray
        }

        # < 8-bit.
        mat <- matrix(seq(0, 255, length.out=1000), 50, 20)
        storage.mode(mat) <- "integer"

        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_UINT8")
        expect_null(out$placeholder)
        out <- alabaster.matrix:::optimize_storage(fun(mat - 128L))
        expect_equal(out$type, "H5T_NATIVE_INT8")
        out <- alabaster.matrix:::optimize_storage(fun(mat + 1L))
        expect_equal(out$type, "H5T_NATIVE_UINT16")

        # < 16-bit.
        mat <- matrix(seq(0, 65535, length.out=1000), 50, 20)
        storage.mode(mat) <- "integer"

        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_UINT16")
        expect_null(out$placeholder)
        out <- alabaster.matrix:::optimize_storage(fun(mat - 32768L))
        expect_equal(out$type, "H5T_NATIVE_INT16")
        out <- alabaster.matrix:::optimize_storage(fun(mat + 1L))
        expect_equal(out$type, "H5T_NATIVE_INT32")

        # < 32-bit.
        mat <- matrix(seq(0, .Machine$integer.max, length.out=1000), 50, 20)
        storage.mode(mat) <- "integer"

        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_INT32")
        expect_null(out$placeholder)
    }
})

test_that("storage optimization works for integers, plus missing values", {
    old <- getAutoBlockSize()
    setAutoBlockSize(400)
    on.exit(setAutoBlockSize(old))

    for (i in 1:2) {
        fun <- identity
        if (i == 1L) {
            fun <- DelayedArray
        }

        # < 8-bit, unsigned.
        mat <- matrix(seq(0, 255, length.out=1000), 50, 20)
        mat[1000] <- NA
        storage.mode(mat) <- "integer" # Note that cast to int does a truncation, so there's only ever one value at a non-zero extreme.

        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_UINT8")
        expect_equal(out$placeholder, 255L)

        mat[999] <- 255L
        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_UINT16")
        expect_equal(out$placeholder, 65535L)

        # 8-bit, signed.
        mat <- matrix(seq(-128, 127, length.out=1000), 50, 20)
        mat[1] <- NA
        storage.mode(mat) <- "integer"

        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_INT8")
        expect_equal(out$placeholder, -128L)

        mat[2] <- -128L
        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_INT16")
        expect_equal(out$placeholder, -32768L)

        # < 16-bit, unsigned.
        mat <- matrix(seq(0, 65535, length.out=1000), 50, 20)
        mat[1000] <- NA
        storage.mode(mat) <- "integer"

        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_UINT16")
        expect_equal(out$placeholder, 65535)

        mat[999] <- 65535L 
        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_INT32")
        expect_equal(out$placeholder, NA_integer_)

        # 16-bit, signed.
        mat <- matrix(seq(-32768, 32767, length.out=1000), 50, 20)
        mat[1] <- NA
        storage.mode(mat) <- "integer"

        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_INT16")
        expect_equal(out$placeholder, -32768L)

        mat[2] <- -32768L
        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_INT32")
        expect_equal(out$placeholder, NA_integer_)

        # Everything else.
        mat <- matrix(seq(0, .Machine$integer.max, length.out=1000), 50, 20)
        storage.mode(mat) <- "integer"
        mat[1] <- NA
        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_INT32")
        expect_equal(out$placeholder, NA_integer_)
    }
})

test_that("storage optimization works for integer-like doubles, no missing values", {
    old <- getAutoBlockSize()
    setAutoBlockSize(400)
    on.exit(setAutoBlockSize(old))

    for (i in 1:2) {
        fun <- identity
        if (i == 1L) {
            fun <- DelayedArray
        }

        # < 8-bit.
        mat <- round(matrix(seq(0, 255, length.out=1000), 50, 20))

        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_UINT8")
        expect_null(out$placeholder)
        out <- alabaster.matrix:::optimize_storage(fun(mat - 128L))
        expect_equal(out$type, "H5T_NATIVE_INT8")
        out <- alabaster.matrix:::optimize_storage(fun(mat + 1L))
        expect_equal(out$type, "H5T_NATIVE_UINT16")

        # < 16-bit.
        mat <- round(matrix(seq(0, 65535, length.out=1000), 50, 20))

        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_UINT16")
        expect_null(out$placeholder)
        out <- alabaster.matrix:::optimize_storage(fun(mat - 32768L))
        expect_equal(out$type, "H5T_NATIVE_INT16")
        out <- alabaster.matrix:::optimize_storage(fun(mat + 1L))
        expect_equal(out$type, "H5T_NATIVE_UINT32")

        # < 32-bit.
        mat <- round(matrix(seq(0, 2^32-1, length.out=1000), 50, 20))

        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_UINT32")
        expect_null(out$placeholder)
        out <- alabaster.matrix:::optimize_storage(fun(mat - 2^31))
        expect_equal(out$type, "H5T_NATIVE_INT32")
        out <- alabaster.matrix:::optimize_storage(fun(mat + 1L))
        expect_equal(out$type, "H5T_NATIVE_DOUBLE")
    }
})

test_that("storage optimization works for integer-like doubles, plus missing values", {
    old <- getAutoBlockSize()
    setAutoBlockSize(400)
    on.exit(setAutoBlockSize(old))

    for (i in 1:2) {
        fun <- identity
        if (i == 1L) {
            fun <- DelayedArray
        }

        # < 8-bit, unsigned
        mat <- round(matrix(seq(0, 255, length.out=1000), 50, 20))
        mat[mat == 255] <- NA

        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_UINT8")
        expect_equal(out$placeholder, 255)

        mat[999] <- 255
        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_UINT16")
        expect_equal(out$placeholder, 65535)

        # < 8-bit, signed
        mat <- round(matrix(seq(-128, 127, length.out=1000), 50, 20))
        mat[mat==-128] <- NA

        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_INT8")
        expect_equal(out$placeholder, -128)

        mat[2] <- -128 
        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_INT16")
        expect_equal(out$placeholder, -32768)

        # < 16-bit, unsigned
        mat <- round(matrix(seq(0, 65535, length.out=1000), 50, 20))
        mat[mat==65535] <- NA

        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_UINT16")
        expect_equal(out$placeholder, 65535)

        mat[999] <- 65535
        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_UINT32")
        expect_equal(out$placeholder, 2^32-1)

        # < 16-bit, signed 
        mat <- round(matrix(seq(-2^15, 2^15-1, length.out=1000), 50, 20))
        mat[mat==-32768] <- NA

        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_INT16")
        expect_equal(out$placeholder, -32768)

        mat[2] <- -32768
        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_INT32")
        expect_equal(out$placeholder, -2^31)

        # < 32-bit, unsigned
        mat <- round(matrix(seq(0, 2^32-1, length.out=1000), 50, 20))
        mat[mat==2^32-1] <- NA

        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_UINT32")
        expect_equal(out$placeholder, 2^32-1)

        mat[999] <- 2^32-1
        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_DOUBLE")
        expect_equal(out$placeholder, NaN)

        # < 32-bit, signed
        mat <- round(matrix(seq(-2^31, 2^31-1, length.out=1000), 50, 20))
        mat[mat == -2^31] <- NA

        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_INT32")
        expect_equal(out$placeholder, -2^31)

        mat[999] <- -2^31
        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_DOUBLE")
        expect_equal(out$placeholder, NaN)
    }
})

test_that("storage optimization works non-integer doubles", {
    old <- getAutoBlockSize()
    setAutoBlockSize(400)
    on.exit(setAutoBlockSize(old))

    for (i in 1:2) {
        fun <- identity
        if (i == 1L) {
            fun <- DelayedArray
        }

        mat <- matrix(rnorm(1000), 50, 20)
        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_DOUBLE")
        expect_null(out$placeholder)

        # Running through the gamut of missing value placeholders.
        mat <- matrix(rnorm(1000), 50, 20)
        mat[1] <- NA
        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_DOUBLE")
        expect_equal(out$placeholder, NaN)

        mat[2] <- NaN
        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_DOUBLE")
        expect_equal(out$placeholder, Inf)

        mat[3] <- Inf
        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_DOUBLE")
        expect_equal(out$placeholder, -Inf)

        mat[4] <- -Inf
        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_DOUBLE")
        expect_equal(out$placeholder, alabaster.matrix:::lowest_double())

        mat[5] <- alabaster.matrix:::lowest_double()
        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_DOUBLE")
        expect_equal(out$placeholder, alabaster.matrix:::highest_double())

        mat[6] <- alabaster.matrix:::highest_double()
        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_DOUBLE")
        expect_true(!any(mat==out$placeholder, na.rm=TRUE))
    }
})

test_that("storage optimization works for strings", {
    old <- getAutoBlockSize()
    setAutoBlockSize(400)
    on.exit(setAutoBlockSize(old))

    for (i in 1:2) {
        fun <- identity
        if (i == 1L) {
            fun <- DelayedArray
        }

        mat <- matrix(sample(LETTERS, 1000, replace=TRUE), 50, 20)
        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(H5Tget_size(out$type), 1)
        expect_equal(H5Tget_cset(out$type), 0) # aka ASCII
        expect_equal(out$placeholder, NULL)

        mat[1] <- NA
        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(H5Tget_size(out$type), 2)
        expect_equal(out$placeholder, "NA")

        mat[2] <- "NA"
        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(H5Tget_size(out$type), 3)
        expect_equal(out$placeholder, "_NA")

        mat[3] <- "_NA"
        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(H5Tget_size(out$type), 4)
        expect_equal(out$placeholder, "__NA")

        # Correct size determination.
        mat <- matrix(sample(LETTERS, 1000, replace=TRUE), 50, 20)
        mat[1] <- "Aaron"
        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(H5Tget_size(out$type), 5) 

        mat <- matrix(sample(LETTERS, 1000, replace=TRUE), 50, 20)
        mat[1000] <- "Aaron"
        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(H5Tget_size(out$type), 5) 

        mat <- matrix("", 50, 20)
        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(H5Tget_size(out$type), 1) 

        # Checking for correct behavior in all-NA cases.
        mat <- matrix(NA_character_, 50, 20)
        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(H5Tget_size(out$type), 2) 
        expect_identical(out$placeholder, "NA")

        # Handles UTF-8. Assumes that this file was sourced as UTF-8.
        mat <- matrix("α", 50, 20)
        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(H5Tget_size(out$type), 2) 
        if (Encoding(mat[1]) == "UTF-8") {
            expect_equal(H5Tget_cset(out$type), 1) # aka UTF-8
        }
    }
})

test_that("storage optimization works for booleans", {
    old <- getAutoBlockSize()
    setAutoBlockSize(400)
    on.exit(setAutoBlockSize(old))

    for (i in 1:2) {
        fun <- identity
        if (i == 1L) {
            fun <- DelayedArray
        }

        mat <- matrix(c(TRUE, FALSE), 50, 20)
        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_INT8")
        expect_null(out$placeholder)

        mat[1] <- NA
        out <- alabaster.matrix:::optimize_storage(fun(mat))
        expect_equal(out$type, "H5T_NATIVE_INT8")
        expect_equal(out$placeholder, -1L)
    }
})

test_that("storage optimization works for sparse objects", {
    old <- getAutoBlockSize()
    setAutoBlockSize(400)
    on.exit(setAutoBlockSize(old))

    for (i in 1:3) {
        if (i == 1L) {
            fun <- function(x) as(x, "sparseMatrix")
        } else if (i == 2L) {
            fun <- function(x) as(x, "SVT_SparseMatrix")
        } else {
            fun <- function(x) DelayedArray(as(x, "COO_SparseArray"))
        }

        # Integer.
        {
             y <- matrix(0L, nrow=100, ncol=10)
             y[1+sample(999, 100)] <- 10000L
             y[1] <- NA

             out <- alabaster.matrix:::optimize_storage(fun(y))
             expect_identical(out$type, "H5T_NATIVE_UINT16")
             expect_equal(out$placeholder, 2^16-1)
             expect_identical(out$size, 101L)
        }

        # Boolean
        {
             y <- matrix(FALSE, nrow=100, ncol=10)
             y[1+sample(999, 50)] <- TRUE 
             y[1] <- NA

             out <- alabaster.matrix:::optimize_storage(fun(y))
             expect_identical(out$type, "H5T_NATIVE_INT8")
             expect_equal(out$placeholder, -1L)
             expect_identical(out$size, 51L)
        }

        # Number
        {
             y <- matrix(0, nrow=100, ncol=10)
             y[1+sample(999, 200)] <- runif(200)
             y[1] <- NA

             out <- alabaster.matrix:::optimize_storage(fun(y))
             expect_identical(out$type, "H5T_NATIVE_DOUBLE")
             expect_equal(out$placeholder, NaN)
             expect_identical(out$size, 201L)
        }
    }
})
