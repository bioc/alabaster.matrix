# library(testthat); library(alabaster.matrix); source("setup.R"); source("test-saveHDF5Array.R")

library(HDF5Array)
test_that("HDF5Arrays are not inadvertently saved as DelayedArrays", {
    mat <- as(array(rpois(1000, 10), c(50, 20)), "HDF5Array")
    dir <- tempfile()
    saveObject(mat, dir)

    obj <- readObjectFile(dir)
    expect_identical(obj$type, "dense_array")

    roundtrip <- readObject(dir)
    expect_identical(as.matrix(roundtrip), as.matrix(mat))
})

test_that("H5SparseMatrices are not inadvertently saved as DelayedArrays", {
    # Easiest to just create a H5SparseMatrix by saving an ordinary matrix.
    mat0 <- Matrix::rsparsematrix(100, 200, 0.1)
    dir0 <- tempfile()
    saveObject(mat0, dir0)
    mat1 <- readObject(dir0) 

    mat <- DelayedArray(mat1@seed@seed)
    expect_s4_class(mat, "H5SparseMatrix")
    dir <- tempfile()
    saveObject(mat, dir)

    obj <- readObjectFile(dir)
    expect_identical(obj$type, "compressed_sparse_matrix")

    roundtrip <- readObject(dir)
    expect_identical(as.matrix(roundtrip), as.matrix(mat))
})
