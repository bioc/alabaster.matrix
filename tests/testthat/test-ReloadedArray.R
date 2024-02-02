# This checks the ReloadedArray, in particular the saveObject method.
# library(testthat); library(alabaster.matrix); source("test-ReloadedArray.R")

arr <- array(rpois(10000, 10), c(50, 20, 10))
dir <- tempfile()
saveObject(arr, dir)
obj <- readObject(dir)

x <- Matrix::rsparsematrix(1000, 200, 0.1)
dir2 <- tempfile()
saveObject(x, dir2)
obj2 <- readObject(dir2)

test_that("ReloadedArrays work correctly", {
    expect_s4_class(obj, "ReloadedArray")
    expect_identical(BiocGenerics::path(obj), dir)
    expect_identical(dim(obj), dim(arr))
    expect_identical(as.array(obj), arr)
    expect_false(is_sparse(obj))

    expect_s4_class(obj2, "ReloadedMatrix")
    expect_identical(BiocGenerics::path(obj2), dir2)
    expect_identical(dim(obj2), dim(x))
    expect_identical(as(obj2, "dgCMatrix"), x)
    expect_true(is_sparse(obj2))
})

test_that("ReloadedArrays save correctly", {
    tmp <- tempfile()
    saveObject(obj, tmp, reloadedarray.reuse.files="none")
    expect_identical(as.array(readObject(tmp)), arr)

    if (.Platform$OS.type=="unix") { 
        # This test just doesn't seem to work on Windows. Either the symlink
        # doesn't form properly, or the symlink path isn't what we expect.
        tmp <- tempfile()
        saveObject(obj, tmp, reloadedarray.reuse.files="symlink")
        expect_identical(as.array(readObject(tmp)), arr)
        expect_identical(Sys.readlink(file.path(tmp, "array.h5")), file.path(dir, "array.h5"))
    }

    # file.info() doesn't report the inode number so we don't have an easy way
    # to distinguish between hard links and a copy. Oh well.
    tmp <- tempfile()
    saveObject(obj, tmp, reloadedarray.reuse.files="link")
    expect_identical(as.array(readObject(tmp)), arr)

    tmp <- tempfile()
    saveObject(obj, tmp, reloadedarray.reuse.files="copy")
    expect_identical(as.array(readObject(tmp)), arr)
})
