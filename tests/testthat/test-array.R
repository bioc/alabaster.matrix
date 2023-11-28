# This tests the saveObject generic for arrays.
# library(testthat); library(alabaster.matrix); source("test-array.R")

library(DelayedArray)
experiment <- "rnaseq"
assay <- "counts"

arr <- array(rpois(10000, 10), c(50, 20, 10))
dimnames(arr) <- list(
   paste0("GENE_", seq_len(nrow(arr))),
   letters[1:20],
   NULL
)

test_that("saveObject works as expected", {
    tmp <- tempfile()
    saveObject(arr, tmp)
    expect_identical(as.array(readObject(tmp)), arr)

    # Works when it's officially integer.
    copy <- arr
    storage.mode(copy) <- "integer"
    tmp <- tempfile()
    saveObject(copy, tmp)
    expect_identical(as.array(readObject(tmp)), copy)

    # Works without dimnames.
    copy <- arr
    dimnames(copy) <- NULL
    tmp <- tempfile()
    saveObject(copy, tmp)
    expect_identical(as.array(readObject(tmp)), copy)
})

test_that("saveObject type optimization works as expected for integers", {
    # Small unsigned integers
    mat <- matrix(sample(255, 1000, replace=TRUE), 40, 25)
    tmp <- tempfile()
    saveObject(mat, tmp)
    expect_identical(as.matrix(readObject(tmp)), mat)

    storage.mode(mat) <- "integer"
    tmp <- tempfile()
    saveObject(mat, tmp)
    expect_identical(as.matrix(readObject(tmp)), mat)

    mat[100] <- NA
    tmp <- tempfile()
    saveObject(mat, tmp)
    expect_identical(as.matrix(readObject(tmp)), mat)

    # Small signed integers
    mat <- matrix(sample(255, 1000, replace=TRUE) - 128, 40, 25)
    tmp <- tempfile()
    saveObject(mat, tmp)
    expect_identical(as.matrix(readObject(tmp)), mat)

    storage.mode(mat) <- "integer"
    tmp <- tempfile()
    saveObject(mat, tmp)
    expect_identical(as.matrix(readObject(tmp)), mat)

    mat[100] <- NA
    tmp <- tempfile()
    saveObject(mat, tmp)
    expect_identical(as.matrix(readObject(tmp)), mat)

    # Large integers
    mat <- trunc(matrix(runif(1000, -1e6, 1e6), 40, 25))
    tmp <- tempfile()
    saveObject(mat, tmp)
    expect_identical(as.matrix(readObject(tmp)), mat)

    storage.mode(mat) <- "integer"
    tmp <- tempfile()
    saveObject(mat, tmp)
    expect_identical(as.matrix(readObject(tmp)), mat)

    mat[100] <- NA
    tmp <- tempfile()
    saveObject(mat, tmp)
    expect_identical(as.matrix(readObject(tmp)), mat)
})

test_that("saveObject works as expected for floats", {
    mat <- matrix(rnorm(1000, -1e6, 1e6), 40, 25)
    tmp <- tempfile()
    saveObject(mat, tmp)
    expect_identical(as.matrix(readObject(tmp)), mat)

    mat[100] <- NA
    tmp <- tempfile()
    saveObject(mat, tmp)
    expect_identical(as.matrix(readObject(tmp)), mat)

    mat[101] <- NaN
    tmp <- tempfile()
    saveObject(mat, tmp)
    expect_identical(as.matrix(readObject(tmp)), mat)

    mat[102] <- Inf
    mat[103] <- -Inf
    tmp <- tempfile()
    saveObject(arr, tmp)
    expect_identical(as.array(readObject(tmp)), arr)
})

test_that("saveObject works as expected for logicals", {
    mat <- matrix(rbinom(1000, 1, 0.5) == 1, 40, 25)
    tmp <- tempfile()
    saveObject(mat, tmp)
    expect_identical(as.matrix(readObject(tmp)), mat)

    mat[100] <- NA
    tmp <- tempfile()
    saveObject(mat, tmp)
    expect_identical(as.matrix(readObject(tmp)), mat)
})

test_that("saveObject works as expected for strings", {
    mat <- matrix(sample(LETTERS, 1000, replace=TRUE), 40, 25)
    tmp <- tempfile()
    saveObject(mat, tmp)
    expect_identical(as.matrix(readObject(tmp)), mat)

    mat[100] <- NA
    tmp <- tempfile()
    saveObject(mat, tmp)
    expect_identical(as.matrix(readObject(tmp)), mat)
})

test_that("saveObject diverts correctly with pristine dense DelayedArrays", {
    x <- DelayedArray(arr)
    expect_true(isPristine(x))

    tmp <- tempfile()
    saveObject(x, tmp)
    expect_identical(as.array(readObject(tmp)), x@seed)
})

test_that("saveObject works correctly with dense block processing", {
    x <- DelayedArray(arr) * 1L
    expect_false(isPristine(x))

    tmp <- tempfile(fileext=".h5")
    local({
        old <- getAutoBlockSize()
        setAutoBlockSize(20000)
        on.exit(setAutoBlockSize(old))
        saveObject(x, tmp)
    })

    roundtrip <- readObject(tmp)
    expect_identical(as.array(roundtrip), arr)
})

test_that("reading dense arrays work with non-default NA placeholders", {
    tmp <- tempfile()
    saveObject(arr, tmp)
    first <- arr[1]

    library(rhdf5)
    local({ 
        fhandle <- H5Fopen(file.path(tmp, "array.h5"))
        on.exit(H5Fclose(fhandle), add=TRUE, after=FALSE)
        ghandle <- H5Gopen(fhandle, "dense_array")
        on.exit(H5Gclose(ghandle), add=TRUE, after=FALSE)
        dhandle <- H5Dopen(ghandle, "data")
        on.exit(H5Dclose(dhandle), add=TRUE, after=FALSE)
        alabaster.base::h5_write_attribute(dhandle, "missing-value-placeholder", first, scalar=TRUE)
    })

    ref <- arr
    ref[ref == first] <- NA
    expect_identical(ref, as.array(readObject(tmp)))
})
