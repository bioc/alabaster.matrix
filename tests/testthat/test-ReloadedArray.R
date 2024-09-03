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
    expect_identical(BiocGenerics::path(obj), normalizePath(dir))
    expect_identical(dim(obj), dim(arr))
    expect_identical(extract_array(obj, vector("list", 3)), arr)
    expect_identical(as.array(obj), arr)
    expect_false(is_sparse(obj))

    expect_s4_class(obj2, "ReloadedMatrix")
    expect_identical(BiocGenerics::path(obj2), normalizePath(dir2))
    expect_identical(dim(obj2), dim(x))
    expect_identical(as(extract_sparse_array(obj2, vector("list", 2)), "dgCMatrix"), x)
    expect_true(is_sparse(obj2))

    # Coercions work as expected.
    expect_identical(as(obj2, "dgCMatrix"), x)
    expect_identical(as(obj2@seed, "dgCMatrix"), x)
})

test_that("ReloadedArrays store the absolute path", {
    cwd <- getwd()
    setwd(dirname(dir))
    on.exit(setwd(cwd))

    obj <- readObject(basename(dir))
    expect_identical(BiocGenerics::path(obj), normalizePath(dir))
    expect_identical(as.array(obj), arr)
})

test_that("ReloadedArrays save correctly", {
    original <- tempfile()
    saveObject(obj, original, ReloadedArray.reuse.files="none")
    expect_identical(as.array(readObject(original)), arr)

    if (.Platform$OS.type=="unix") { 
        # This test just doesn't seem to work on Windows. Either the symlink
        # doesn't form properly, or the symlink path isn't what we expect.
        tmp <- tempfile()
        saveObject(obj, tmp, ReloadedArray.reuse.files="symlink")
        expect_identical(as.array(readObject(tmp)), arr)
        link.dest <- Sys.readlink(file.path(tmp, "array.h5"))
        expect_true(startsWith(link.dest, "/"))
        expect_identical(normalizePath(link.dest), normalizePath(file.path(dir, "array.h5")))
        expect_identical(as.array(readObject(tmp)), arr)

        # Relative symlinks also work as expected.
        tmp <- tempfile()
        saveObject(obj, tmp, ReloadedArray.reuse.files="relsymlink")
        expect_identical(as.array(readObject(tmp)), arr)
        link.dest <- Sys.readlink(file.path(tmp, "array.h5"))
        expect_true(startsWith(link.dest, ".."))
        expect_identical(normalizePath(file.path(tmp, link.dest)), normalizePath(file.path(dir, "array.h5")))
        expect_identical(as.array(readObject(tmp)), arr)

        # Trying in a more deeply nested target directory.
        tmp0 <- tempfile()
        dir.create(tmp0)
        tmp <- tempfile(tmpdir=tmp0)
        saveObject(obj, tmp, ReloadedArray.reuse.files="relsymlink")
        expect_identical(as.array(readObject(tmp)), arr)
        link.dest <- Sys.readlink(file.path(tmp, "array.h5"))
        expect_true(startsWith(link.dest, "../.."))
        expect_identical(normalizePath(file.path(tmp, link.dest)), normalizePath(file.path(dir, "array.h5")))
        expect_identical(as.array(readObject(tmp)), arr)
    }

    # file.info() doesn't report the inode number so we don't have an easy way
    # to distinguish between hard links and a copy. Oh well.
    tmp <- tempfile()
    saveObject(obj, tmp, ReloadedArray.reuse.files="link")
    expect_identical(as.array(readObject(tmp)), arr)

    tmp <- tempfile()
    saveObject(obj, tmp, ReloadedArray.reuse.files="copy")
    expect_identical(as.array(readObject(tmp)), arr)
})
