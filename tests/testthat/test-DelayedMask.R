# Tests that the delayed masking works as expected.
# library(testthat); library(alabaster.matrix); source("test-DelayedMask.R")

library(DelayedArray)
set.seed(100000)

test_that("basic tests for DelayedMask", {
    original <- matrix(rpois(400, lambda=2), ncol=10)
    storage.mode(original) <- "integer"
    masked <- DelayedArray(DelayedMask(original, 0))
    expect_identical(dim(original), dim(masked))
    expect_identical(type(original), type(masked))
    expect_false(is_sparse(masked))

    ref <- original
    ref[ref == 0] <- NA
    expect_identical(ref, as.matrix(masked))

    # Other methods work as expected.
    y <- original
    rownames(y) <- sprintf("GENE_%i", seq_len(nrow(y)))
    masked <- DelayedArray(DelayedMask(y, 0))
    expect_identical(rownames(masked), rownames(y))
})

test_that("sparse tests for DelayedMask", {
    original <- round(Matrix::rsparsematrix(100, 20, density=0.2) * 10)
    masked <- DelayedArray(DelayedMask(original, 1))
    expect_true(is_sparse(masked))

    ref <- as.matrix(original)
    dimnames(ref) <- list(NULL, NULL)
    ref[ref == 1] <- NA
    expect_identical(ref, as.matrix(masked))
    dimnames(ref) <- NULL

#    spmat <- extract_sparse_array(masked, list(NULL, NULL))
#    expect_identical(ref, as.matrix(spmat))

    spseed <- OLD_extract_sparse_array(masked, list(NULL, NULL))
    expect_identical(ref, as.matrix(spseed))

    # Behaves correctly if zero is the placeholder.
    masked <- DelayedArray(DelayedMask(original, 0))
    expect_false(is_sparse(masked))

    ref <- as.matrix(original)
    ref[ref == 0] <- NA
    dimnames(ref) <- list(NULL, NULL)
    expect_identical(ref, as.matrix(masked))
})

test_that("DelayedMask works with special values", {
    original <- matrix(rpois(400, lambda=2), ncol=10)
    expect_identical(DelayedMask(original, NA), original)

    masked <- DelayedMask(original, NA, force=TRUE)
    expect_identical(original, as.matrix(masked))

    # Also works with NaN values, whether the placeholder is NaN or not.
    original[1,1] <- NaN

    masked <- DelayedMask(original, 1)
    ref <- as.matrix(original)
    ref[ref == 1] <- NA
    expect_identical(ref, as.matrix(masked))

    masked <- DelayedMask(original, NaN)
    ref <- as.matrix(original)
    ref[is.nan(ref)] <- NA
    expect_identical(ref, as.matrix(masked))

    # Same with a different type.
    copy <- original
    storage.mode(copy) <- "integer"
    expect_identical(DelayedMask(copy, NA), copy)
    expect_identical(DelayedMask(copy, NaN), copy)
})
