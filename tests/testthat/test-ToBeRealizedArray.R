# library(testthat); library(alabaster.matrix); source("test-ToBeRealizedArray.R")

arr <- array(rpois(10000, 10), c(50, 20, 10))
x <- Matrix::rsparsematrix(1000, 200, 0.1)

test_that("ToBeRealizedArrays work correctly", {
    obj <- ToBeRealizedArray(arr)
    expect_s4_class(obj, "ToBeRealizedArray")
    expect_identical(dim(obj), dim(arr))
    expect_identical(extract_array(obj, vector("list", 3)), arr)
    expect_identical(as.array(obj), arr)
    expect_false(is_sparse(obj))

    obj2 <- ToBeRealizedArray(x)
    expect_s4_class(obj2, "ToBeRealizedArray")
    expect_identical(dim(obj2), dim(x))
    expect_identical(as(extract_sparse_array(obj2, vector("list", 2)), "dgCMatrix"), x)
    expect_true(is_sparse(obj2))

    # Coercions work as expected.
    expect_identical(as(obj2, "dgCMatrix"), x)
    expect_identical(as(obj2@seed, "dgCMatrix"), x)
})

test_that("ToBeRealizedArrays save correctly", {
    foo <- log1p(DelayedArray(arr))
    bar <- (DelayedArray(x) * 10) ^2

    tmp <- tempfile()
    saveObject(
        list(
            foo1=foo,
            foo2=ToBeRealizedArray(foo),
            bar1=bar,
            bar2=ToBeRealizedArray(bar)
        ),
        tmp,
        DelayedArray.preserve.ops=TRUE
    )

    expect_identical(readObjectFile(file.path(tmp, "other_contents", "0"))$type, "delayed_array")
    expect_identical(readObjectFile(file.path(tmp, "other_contents", "1"))$type, "dense_array")
    expect_identical(readObjectFile(file.path(tmp, "other_contents", "2"))$type, "delayed_array")
    expect_identical(readObjectFile(file.path(tmp, "other_contents", "3"))$type, "compressed_sparse_matrix")

    roundtrip <- readObject(tmp)
    reffoo <- log1p(arr)
    expect_identical(as.array(roundtrip$foo1), reffoo)
    expect_identical(as.array(roundtrip$foo2), reffoo)
    refbar <- (x * 10) ^ 2
    expect_identical(as(roundtrip$bar1, "dgCMatrix"), refbar)
    expect_identical(as(roundtrip$bar2, "dgCMatrix"), refbar)
})

test_that("ToBeRealizedArrays respect deduplication", {
    foo <- log1p(DelayedArray(arr))
    bar <- (DelayedArray(x) * 10) ^2

    payload <- list(
        foo1 = ToBeRealizedArray(foo),
        nested = list(foo2 = ToBeRealizedArray(foo), bar1 = ToBeRealizedArray(bar)),
        bar2 = ToBeRealizedArray(bar)
    )

    tmp <- tempfile()
    saveObject(
        payload,
        tmp,
        DelayedArray.preserve.ops = TRUE,
        array.dedup.session = createDedupSession(),
        array.dedup.action = ifelse(.Platform$OS.type == "unix", "relsymlink", "copy")
    )

    if (.Platform$OS.type=="unix") { 
        link.dest <- file.path(tmp, "other_contents", "0", "array.h5")
        expect_identical(Sys.readlink(link.dest), "")
        link.dest <- file.path(tmp, "other_contents", "1", "array.h5")
        expect_true(startsWith(Sys.readlink(link.dest), ".."))
        link.dest <- file.path(tmp, "other_contents", "2", "matrix.h5")
        expect_identical(Sys.readlink(link.dest), "")
        link.dest <- file.path(tmp, "other_contents", "3", "matrix.h5")
        expect_true(startsWith(Sys.readlink(link.dest), ".."))
    }

    roundtrip <- readObject(tmp)
    expect_identical(as.array(roundtrip$foo1), as.array(foo))
    expect_identical(as.array(roundtrip$nested$foo2), as.array(foo))
    expect_identical(as(roundtrip$nested$bar1, "dgCMatrix"), as(bar, "dgCMatrix"))
    expect_identical(as(roundtrip$bar2, "dgCMatrix"), as(bar, "dgCMatrix"))
})

test_that("We correctly save ToBeRealizedArraySeeds inside DelayedArrays", {
    X <- ToBeRealizedArray(log1p(DelayedArray(arr))) * 10
    tmp <- tempfile()
    saveObject(X, tmp, DelayedArray.preserve.ops = TRUE)

    roundtrip <- readObject(tmp)
    expect_identical(as.array(X), as.array(roundtrip))

    # The log1p operation is correctly realized.
    internal <- readObject(file.path(tmp, "seeds", "0"))
    expect_identical(log1p(arr), as.array(internal))
})
