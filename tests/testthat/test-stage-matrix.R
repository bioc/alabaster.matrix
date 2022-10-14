# This tests the stageObject generic.
# library(testthat); library(alabaster.matrix); source("test-stage-matrix.R")

experiment <- "rnaseq"
assay <- "counts"

mat <- matrix(rpois(10000, 10), ncol=10)
rownames(mat) <- paste0("GENE_", seq_len(nrow(mat)))
colnames(mat) <- letters[1:10]

test_that("stageObject works as expected for the default method", {
    dir <- tempfile()
    odir <- file.path(dir, experiment)
    dir.create(odir, recursive=TRUE)

    info <- stageObject(mat, dir, file.path(experiment, assay))
    expect_match(info$`$schema`, "hdf5_dense_array")

    # Metadata is valid.
    expect_error(alabaster.base::.writeMetadata(info, dir), NA)

    # Round tripping.
    mat2 <- loadArray(info, project=dir)
    expect_equal(sum(mat2), sum(mat))
    expect_identical(rownames(mat), rownames(mat2))
    expect_identical(colnames(mat), colnames(mat2))
})

test_that("stageObject works as expected without dimnames", {
    dir <- tempfile()
    odir <- file.path(dir, experiment)
    dir.create(odir, recursive=TRUE)

    dimnames(mat) <- NULL
    info <- stageObject(mat, dir, file.path(experiment, assay))

    # Metadata is valid.
    expect_error(alabaster.base::.writeMetadata(info, dir), NA)

    # Round-tripping.
    mat2 <- loadArray(info, project=dir)
    expect_equal(sum(mat2), sum(mat))
    expect_null(rownames(mat2))
    expect_null(colnames(mat2))
})

test_that("stageObject works with sparse matrices", {
    dir <- tempfile()
    odir <- file.path(dir, experiment)
    dir.create(odir, recursive=TRUE)

    mat <- as(matrix(rpois(10000, 1), ncol=10), "dgCMatrix")
    rownames(mat) <- paste0("GENE_", seq_len(nrow(mat)))
    colnames(mat) <- letters[1:10]

    info <- stageObject(mat, dir, file.path(experiment, assay))
    expect_match(info$`$schema`, "sparse_matrix")

    # Metadata is valid.
    expect_error(alabaster.base::.writeMetadata(info, dir), NA)

    # Round tripping.
    mat2 <- loadArray(info, project=dir)
    expect_equal(sum(mat2), sum(mat))
    expect_identical(rownames(mat), rownames(mat2))
    expect_identical(colnames(mat), colnames(mat2))

    # Without dimnames.
    dir <- tempfile()
    odir <- file.path(dir, experiment)
    dir.create(odir, recursive=TRUE)

    mat <- as(matrix(rpois(10000, 1), ncol=10), "dgCMatrix")
    info <- stageObject(mat, dir, file.path(experiment, assay))
    expect_match(info$`$schema`, "sparse_matrix")
    expect_error(alabaster.base::.writeMetadata(info, dir), NA)

    mat2 <- loadArray(info, project=dir)
    expect_null(rownames(mat2))
    expect_null(colnames(mat2))
})

test_that("stageObject works with character matrices and missing values", {
    dir <- tempfile()
    odir <- file.path(dir, experiment)
    dir.create(odir, recursive=TRUE)

    mat <- matrix(sample(LETTERS, 100, replace=TRUE), 25, 4)
    mat[1] <- "NA"
    mat[100] <- NA_character_

    info <- stageObject(mat, dir, file.path(experiment, assay))
    expect_match(info$`$schema`, "hdf5_dense")

    # Metadata is valid.
    expect_error(alabaster.base::.writeMetadata(info, dir), NA)

    # Round-tripping.
    mat2 <- loadArray(info, project=dir)
    expect_identical(mat, as.matrix(mat2))
})

test_that("stageObject works with DelayedMatrices (naive)", {
    dir <- tempfile()
    odir <- file.path(dir, experiment)
    dir.create(odir, recursive=TRUE)

    # Works in the pristine case, where it just extracts the seed.
    library(DelayedArray)
    mat <- as(matrix(rpois(10000, 1), ncol=10), "dgCMatrix")
    x <- DelayedArray(mat)

    info <- stageObject(x, dir, file.path(experiment, assay))
    expect_match(info$`$schema`, "sparse_matrix")
    expect_error(alabaster.base::.writeMetadata(info, dir), NA)

    mat2 <- loadArray(info, project=dir)
    expect_equal(unname(as.matrix(mat2)), unname(as.matrix(x)))
    expect_null(rownames(mat2))
    expect_null(colnames(mat2))

    # Works with assigned dimnames.
    y <- x
    dimnames(y) <- list(seq_len(nrow(y)), LETTERS[1:10])

    info <- stageObject(y, dir, file.path(experiment, "other"))
    expect_error(alabaster.base::.writeMetadata(info, dir), NA)

    mat2 <- loadArray(info, project=dir)
    expect_equal(as.matrix(mat2), as.matrix(y))
})
