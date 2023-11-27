# This tests the stageObject generic for arrays.
# library(testthat); library(alabaster.matrix); source("test-save-array.R")

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

    # Trying in the new world.
    tmp <- tempfile()
    saveObject(arr, tmp)
    expect_identical(as.array(readArray(tmp)), arr)
})

test_that("stageObject works as expected for NA values", {
    dir <- tempfile()
    odir <- file.path(dir, experiment)
    dir.create(odir, recursive=TRUE)

    storage.mode(arr) <- "integer"
    {
        info <- stageObject(arr, dir, file.path(experiment, "assay-0"))
        fpath <- file.path(odir, "assay-0/array.h5")
        out <- rhdf5::h5readAttributes(fpath, "data")
        expect_false("missing-value-placeholder" %in% names(out))

        arr2 <- loadArray(info, project=dir)
        expect_identical(DelayedArray::type(arr2), "integer")
        expect_equal(arr, as.array(arr2))

        tmp <- tempfile()
        saveObject(arr, tmp)
        expect_identical(as.array(readArray(tmp)), arr)
    }

    arr[1] <- NA
    {
        info <- stageObject(arr, dir, file.path(experiment, "assay-1"))
        fpath <- file.path(odir, "assay-1/array.h5")
        out <- rhdf5::h5readAttributes(fpath, "data")
        expect_identical(out[["missing-value-placeholder"]], NA_integer_)

        arr2 <- loadArray(info, project=dir)
        expect_identical(DelayedArray::type(arr2), "integer")
        expect_equal(arr, as.array(arr2))

        tmp <- tempfile()
        saveObject(arr, tmp)
        expect_identical(as.array(readArray(tmp)), arr)
    }

    storage.mode(arr) <- "double"
    arr[2] <- NaN
    arr[3] <- Inf
    arr[4] <- -Inf
    {
        info <- stageObject(arr, dir, file.path(experiment, "assay-2"))
        fpath <- file.path(odir, "assay-2/array.h5")
        out <- rhdf5::h5readAttributes(fpath, "data")
        #expect_identical(out[["missing-value-placeholder"]], alabaster.matrix:::lowest_double())

        arr2 <- loadArray(info, project=dir)
        expect_identical(DelayedArray::type(arr2), "double")
        expect_equal(arr, as.array(arr2))

        tmp <- tempfile()
        saveObject(arr, tmp)
        expect_identical(as.array(readArray(tmp)), arr)
    }
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

    tmp <- tempfile()
    saveObject(arr, tmp)
    expect_identical(as.array(readArray(tmp)), arr)
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

test_that("stageObject works with DelayedArrays with local references", {
    dir <- tempfile()
    odir <- file.path(dir, experiment)
    dir.create(odir, recursive=TRUE)

    info <- stageObject(arr, dir, paste0(experiment, "/raw"))
    alabaster.base::.writeMetadata(info, dir=dir)
    raw.path <- info$path

    # Adding delayed operations.
    library(DelayedArray)
    x <- DelayedArray(arr)
    x <- log2(x / runif(nrow(x)) + 1)

    old <- preserveDelayedOperations(TRUE)
    info <- stageObject(x, dir, file.path(experiment, assay))
    preserveDelayedOperations(old)

    expect_match(info$`$schema`, "hdf5_delayed_array")
    expect_error(alabaster.base::.writeMetadata(info, dir=dir), NA)

    # Manually replacing the seed with a local reference.
    fpath <- file.path(dir, info$path)
    rhdf5::h5delete(fpath, "data/seed/seed/seed")
    rhdf5::h5createGroup(fpath, "data/seed/seed/seed")
    rhdf5::h5write(dim(arr), fpath, "data/seed/seed/seed/dimensions")
    rhdf5::h5write("integer", fpath, "data/seed/seed/seed/type")
    rhdf5::h5write(raw.path, fpath, "data/seed/seed/seed/path")

    {
        fhandle = rhdf5::H5Fopen(fpath, "H5F_ACC_RDWR")
        ghandle = rhdf5::H5Gopen(fhandle, "data/seed/seed/seed")
        rhdf5::h5writeAttribute("array", ghandle, "delayed_type", asScalar=TRUE)
        rhdf5::h5writeAttribute("custom alabaster local array", ghandle, "delayed_array", asScalar=TRUE)
        rhdf5::H5Gclose(ghandle)
        rhdf5::H5Fclose(fhandle)
    }

    # Now seeing if we can load it.
    arr2 <- loadArray(info, project=dir)
    expect_equal(unname(as.array(arr2)), unname(as.array(x)))
})

test_that("stageObject recycles existing HDF5Arrays", {
    dir <- tempfile()
    dir.create(dir, recursive=TRUE)

    library(HDF5Array)
    x <- as(matrix(runif(1000), 10, 100), "HDF5Array")
    rhdf5::h5write("foo", path(x), "bar") # injecting a little something to distinguish it from a rewritten file.

    {
        info <- stageObject(x, dir, "dense")
        expect_false(identical(file.size(path(x)), file.size(file.path(dir, info$path)))) # negative control that triggers the fallback.
        arr2 <- loadArray(info, project=dir)
        expect_equal(unname(as.array(arr2)), unname(as.array(x)))
    }

    old <- recycleHdf5Files(TRUE)
    on.exit(recycleHdf5Files(old))

    {
        info <- stageObject(x, dir, "dense2")
        expect_match(info$`$schema`, "hdf5_dense_array")
        expect_error(alabaster.base::.writeMetadata(info, dir=dir), NA)
        expect_identical(file.size(path(x)), file.size(file.path(dir, info$path))) # same file.

        arr2 <- loadArray(info, project=dir)
        expect_equal(unname(as.array(arr2)), unname(as.array(x)))
    }

    {
        y <- Matrix::rsparsematrix(20, 50, density=0.1)
        tmp <- tempfile(fileext=".h5")
        rhdf5::h5createFile(tmp)
        rhdf5::h5createGroup(tmp, "foo")
        rhdf5::h5write(y@p, tmp, "foo/indptr")
        rhdf5::h5write(y@i, tmp, "foo/indices")
        rhdf5::h5write(y@x, tmp, "foo/data")
        rhdf5::h5write(dim(y), tmp, "foo/shape")
        x <- H5SparseMatrix(tmp, "foo")

        info <- stageObject(x, dir, "sparse")
        expect_match(info$`$schema`, "hdf5_sparse_matrix")
        expect_error(alabaster.base::.writeMetadata(info, dir=dir), NA)
        expect_identical(file.size(path(x)), file.size(file.path(dir, info$path)))

        arr2 <- loadArray(info, project=dir)
        expect_equal(unname(as.array(arr2)), unname(as.array(x)))
    }
})

test_that("reading arrays work with non-default NA placeholders", {
    dir <- tempfile()
    dir.create(dir, recursive=TRUE)

    # Dense case.
    {
        x <- matrix(rpois(1000, 2), 100, 10)
        storage.mode(x) <- "integer"

        info <- stageObject(x, dir, "dense")
        alabaster.base::addMissingPlaceholderAttributeForHdf5(
            file.path(dir, info$path), 
            info$hdf5_dense_array$dataset, 
            0L
        )

        arr2 <- loadArray(info, project=dir)
        ref <- as.matrix(x)
        ref[ref==0L] <- NA
        expect_identical(ref, as.matrix(arr2))

        # Trying in the new world.
        tmp <- tempfile()
        saveObject(x, tmp)
        local({ 
            fhandle <- H5Fopen(file.path(tmp, "array.h5"))
            on.exit(H5Fclose(fhandle), add=TRUE, after=FALSE)
            ghandle <- H5Gopen(fhandle, "dense_array")
            on.exit(H5Gclose(ghandle), add=TRUE, after=FALSE)
            dhandle <- H5Dopen(ghandle, "data")
            on.exit(H5Dclose(dhandle), add=TRUE, after=FALSE)
            alabaster.base::h5_write_attribute(dhandle, "missing-value-placeholder", 0L, scalar=TRUE)
        })
        expect_identical(ref, as.matrix(readArray(tmp)))
    }

    # Sparse case.
    {
        x <- rsparsematrix(100, 20, 0.5)
        first <- x@x[1]

        info <- stageObject(x, dir, "sparse")
        alabaster.base::addMissingPlaceholderAttributeForHdf5(
            file.path(dir, info$path), 
            paste0(info$hdf5_sparse_matrix$group, "/data"), 
            first
        )

        arr2 <- loadArray(info, project=dir)
        ref <- as.matrix(x)
        ref[ref==first] <- NA
        expect_identical(ref, as.matrix(arr2))
    }
})
