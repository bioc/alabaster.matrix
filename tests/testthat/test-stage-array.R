# This tests the stageObject generic for arrays.
# library(testthat); library(alabaster.matrix); source("test-stage-array.R")

experiment <- "rnaseq"
assay <- "counts"

arr <- array(rpois(10000, 10), c(50, 20, 10))
dimnames(arr) <- list(
   paste0("GENE_", seq_len(nrow(arr))),
   letters[1:20],
   NULL
)

test_that("stageObject works as expected for the default method", {
    dir <- tempfile()
    odir <- file.path(dir, experiment)
    dir.create(odir, recursive=TRUE)

    info <- stageObject(arr, dir, file.path(experiment, assay))
    expect_match(info$`$schema`, "hdf5_dense_array")

    # Round tripping.
    arr2 <- loadArray(info, project=dir)
    expect_equal(sum(arr2), sum(arr))
    expect_identical(rownames(arr), rownames(arr2))
    expect_identical(colnames(arr), colnames(arr2))

    # Checking that metadata save works.
    expect_error(alabaster.base::.writeMetadata(info, dir=dir), NA)
})

test_that("stageObject works as expected without dimnames", {
    dir <- tempfile()
    odir <- file.path(dir, experiment)
    dir.create(odir, recursive=TRUE)

    dimnames(arr) <- NULL
    info <- stageObject(arr, dir, file.path(experiment, assay))

    # Round-tripping.
    arr2 <- loadArray(info, project=dir)
    expect_equal(sum(arr2), sum(arr))
    expect_null(rownames(arr2))
    expect_null(colnames(arr2))
})

test_that("stageObject works with DelayedArrays", {
    dir <- tempfile()
    odir <- file.path(dir, experiment)
    dir.create(odir, recursive=TRUE)

    # Works in the pristine case.
    library(DelayedArray)
    dimnames(arr) <- NULL
    x <- DelayedArray(arr)

    info <- stageObject(x, dir, file.path(experiment, assay))
    expect_match(info$`$schema`, "hdf5_dense_array")

    arr2 <- loadArray(info, project=dir)
    expect_equal(unname(as.array(arr2)), unname(as.array(x)))
    expect_null(rownames(arr2))
    expect_null(colnames(arr2))

    # Works with assigned dimnames.
    y <- x
    dimnames(y) <- list(as.character(seq_len(50)), NULL, LETTERS[1:10])

    info <- stageObject(y, dir, file.path(experiment, "other"))
    arr2 <- loadArray(info, project=dir)
    expect_equal(as.array(arr2), as.array(y))
})

test_that("stageObject works with DelayedArrays with preserved operations", {
    dir <- tempfile()
    odir <- file.path(dir, experiment)
    dir.create(odir, recursive=TRUE)

    # Works in the pristine case.
    library(DelayedArray)
    x <- DelayedArray(arr)
    x <- log2(x / runif(nrow(x)) + 1)

    old <- preserveDelayedOperations(TRUE)
    info <- stageObject(x, dir, file.path(experiment, assay))
    preserveDelayedOperations(old)
    expect_match(info$`$schema`, "hdf5_delayed_array")
    expect_error(alabaster.base::.writeMetadata(info, dir=dir), NA)

    arr2 <- loadArray(info, project=dir)
    expect_equal(unname(as.array(arr2)), unname(as.array(x)))
})
