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
    masked <- DelayedMask(original, NA)
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

    # If the placeholder is NA, both NaNs and NAs are considered to be missing,
    # as we don't rely on the payload (i.e., same logic as is.na()).
    masked <- DelayedMask(original, NA)
    ref <- as.matrix(original)
    ref[is.nan(ref)] <- NA
    expect_identical(ref, as.matrix(masked))

    # No problem if the dataset is integer.
    as.int <- original
    storage.mode(as.int) <- "integer"
    masked <- DelayedMask(as.int, NA_integer_)
    expect_identical(as.int, as.matrix(masked))
})

test_that("DelayedMask fails for integers with non-NA placeholders", {
    original <- matrix(NA_integer_, 10, 10)
    masked <- DelayedMask(original, NA)
    expect_identical(original, as.matrix(masked))

    masked <- DelayedMask(original, 10L)
    expect_error(as.matrix(masked), "not yet supported")
})
