# Tests the AmalgamatedArray class.
# library(testthat); library(alabaster.matrix); source("test-AmalgamatedArray.R")

test_that("Amalgamated basics work as expected", {
    first <- Matrix::rsparsematrix(10, 10, 0.1)
    second <- Matrix::rsparsematrix(10, 20, 0.1)

    mat <- AmalgamatedArray(list(foo = first, bar = second), along=2)
    expect_s4_class(mat, "AmalgamatedArray")
    expect_s4_class(mat, "AmalgamatedMatrix")
    expect_identical(as.array(mat), as.array(BiocGenerics::cbind(first, second)))

    expect_identical(componentNames(mat), c("foo", "bar"))
    out <- extractComponents(mat)
    expect_identical(out$foo, first)
    expect_identical(out$bar, second)

    # Same for the other dimension.
    tsecond <- Matrix::t(second)
    mat <- AmalgamatedArray(list(stuff = first, whee = tsecond), along=1)
    expect_identical(as.array(mat), as.array(BiocGenerics::rbind(first, tsecond)))
})

library(alabaster.base)
test_that("Amalgamated staging and loading works as expected", {
    first <- Matrix::rsparsematrix(10, 10, 0.1)
    second <- Matrix::rsparsematrix(10, 20, 0.1)
    mat <- AmalgamatedArray(list(foo = first, bar = second), along=2)

    temp <- tempfile()
    dir.create(temp)
    meta <- stageObject(mat, temp, "amalgam")
    expect_identical(meta[["$schema"]], "amalgamated_array/v1.json")
    expect_identical(meta$amalgamated_array$along, 1L)
    info <- .writeMetadata(meta, temp)

    remeta <- acquireMetadata(temp, info$path)
    roundtrip <- loadObject(remeta, temp)
    expect_s4_class(roundtrip, "AmalgamatedArray")
    expect_identical(as.array(roundtrip), as.array(mat))
})
