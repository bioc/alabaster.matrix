# This tests the behavior of the various HDF4-based seeds.
# library(testthat); library(alabaster.matrix); source("test-storeDelayedObject.R")

library(DelayedArray)
library(rhdf5)

saveDelayed <- function(x, ...) {
    tmp <- tempfile()
    saveObject(x, tmp, delayedarray.dispatch.pristine=FALSE, delayedarray.preserve.ops=TRUE, delayedarray.store.args=list(...))
    stopifnot(readObjectFile(tmp)$type == "delayed_array")
    tmp
}

loadDelayed <- function(path, ...) {
    raw <- readObject(path, delayedarray.reload.args=list(...))
    DelayedArray(raw@seed@seed)
}

#######################################################
#######################################################

test_that("ConstantArrays are saved correctly", {
    # Behaves correctly for logical values.
    X <- ConstantArray(dim=c(12, 7))
    temp <- saveDelayed(X)

    out <- loadDelayed(temp)
    expect_identical(X, out)

    # Behaves correctly for string values.
    X <- ConstantArray(dim=c(12, 7), value="foo")
    temp <- saveDelayed(X)

    out <- loadDelayed(temp)
    expect_identical(X, out)

    # Behaves correctly for numeric values.
    X <- ConstantArray(dim=c(12, 7), value=2L)
    temp <- saveDelayed(X)

    out <- loadDelayed(temp)
    expect_identical(X, out)
})

test_that("ConstantArrays are still saved correctly after some deep nesting", {
    X <- ConstantArray(dim=c(12, 7), value=23)
    Y <- DelayedArray(matrix(runif(60), nrow=12))
    Z <- cbind(X, Y)

    temp <- saveDelayed(Z)
    out <- loadDelayed(temp)
    expect_identical(Z, out)
})

test_that("ConstantArrays behave correctly with NAs", {
    X <- ConstantArray(dim=c(12, 7), value=NA)
    temp <- saveDelayed(X)
    out <- loadDelayed(temp)
    expect_identical(X, out)

    # Trying a non-Default NA.
    X <- ConstantArray(dim=c(12, 7), value=1.2)
    temp <- saveDelayed(X)

    library(rhdf5)
    (function() {
        fhandle <- H5Fopen(file.path(temp, "array.h5"))
        on.exit(H5Fclose(fhandle), add=TRUE)
        ghandle <- H5Gopen(fhandle, "delayed_array")
        on.exit(H5Gclose(ghandle), add=TRUE)
        dhandle <- H5Dopen(ghandle, "value")
        on.exit(H5Dclose(dhandle), add=TRUE)
        h5_write_attribute(dhandle, "missing_placeholder", 1.2, type="H5T_NATIVE_DOUBLE", scalar=TRUE)
    })()

    out <- loadDelayed(temp)
    expect_identical(ConstantArray(c(12, 7), value=NA_real_), out)
})

#######################################################
#######################################################

test_that("dense arrays are saved correctly", {
    X <- DelayedArray(matrix(rpois(60, 5), ncol=20))
    temp <- saveDelayed(X)
    roundtrip <- loadDelayed(temp)
    expect_identical(X, roundtrip)

    X <- DelayedArray(Matrix::Matrix(runif(60), ncol=20))
    temp <- saveDelayed(X)
    roundtrip <- loadDelayed(temp)
    expect_identical(as.array(X), as.array(roundtrip@seed))

    X <- DelayedArray(Matrix::Matrix(runif(60) > 0.5, ncol=20, sparse=FALSE))
    temp <- saveDelayed(X)
    roundtrip <- loadDelayed(temp)
    expect_identical(as.array(X), as.array(roundtrip@seed))

    # Check that this works with other dense types.
    X <- DelayedArray(Matrix::triu(Matrix::Matrix(runif(100), ncol=10)))
    expect_s4_class(X@seed, "dtrMatrix")
    temp <- saveDelayed(X)
    roundtrip <- loadDelayed(temp)
    expect_identical(as.array(X), roundtrip@seed)
})

test_that("dense array saving handles dimnames", {
    X <- matrix(rpois(60, 5), ncol=20)
    rownames(X) <- LETTERS[1:3]
    colnames(X) <- letters[1:20]

    temp <- saveDelayed(DelayedArray(X))
    roundtrip <- loadDelayed(temp)
    expect_identical(X, roundtrip@seed)
})

test_that("dense array saving handles a bit of transposition", {
    X <- DelayedArray(matrix(rpois(60, 5), ncol=20))
    temp <- saveDelayed(X)

    (function() {
        fhandle <- H5Fopen(file.path(temp, "array.h5"))
        on.exit(H5Fclose(fhandle), add=TRUE, after=FALSE)
        ghandle <- H5Gopen(fhandle, "delayed_array")
        on.exit(H5Gclose(ghandle), add=TRUE, after=FALSE)
        H5Ldelete(ghandle, "native")
        h5_write_vector(ghandle, "native", 1, type="H5T_NATIVE_UINT8", scalar=TRUE)
    })()

    roundtrip <- loadDelayed(temp)
    expect_identical(t(X@seed), roundtrip@seed)
})

test_that("dense arrays are saved to external arrays if requested", {
    X <- DelayedArray(matrix(rpois(60, 5), ncol=20))
    temp <- saveDelayed(X, save.external.array=TRUE)
    expect_true(file.exists(file.path(temp, "seeds", 0)))
    roundtrip <- loadDelayed(temp, custom.takane.realize=TRUE)
    expect_identical(X, roundtrip)

    X <- DelayedArray(Matrix::Matrix(runif(60), ncol=20))
    temp <- saveDelayed(X, save.external.array=TRUE)
    expect_true(file.exists(file.path(temp, "seeds", 0)))
    roundtrip <- loadDelayed(temp, custom.takane.realize=TRUE)
    expect_identical(unname(as.array(X)), as.array(roundtrip@seed))
})

#######################################################
#######################################################

test_that("sparse matrices are saved correctly", {
    X <- DelayedArray(Matrix::rsparsematrix(20, 30, 0.2))
    temp <- saveDelayed(X)
    roundtrip <- loadDelayed(temp)
    expect_s4_class(roundtrip@seed, "SVT_SparseMatrix")
    expect_identical(unname(as.matrix(X)), as.matrix(roundtrip))

    # Checking that it works in CSR mode. 
    X <- Matrix::rsparsematrix(20, 30, 0.2)
    X <- as(X > 0, "RsparseMatrix")
    Z <- DelayedArray(X)
    temp <- saveDelayed(Z)
    roundtrip <- loadDelayed(temp)
    expect_s4_class(roundtrip@seed, "SVT_SparseMatrix")
    expect_identical(unname(as.matrix(X)), as.matrix(roundtrip))

    # Works with SVT sparse matrix inputs.
    X <- matrix(rpois(100, 0.5) * 10L, ncol=25)
    X <- as(X, "SVT_SparseMatrix")
    Z <- DelayedArray(X)
    temp <- saveDelayed(Z)
    roundtrip <- loadDelayed(temp)
    expect_s4_class(roundtrip@seed, "SVT_SparseMatrix")
    expect_identical(unname(as.matrix(X)), as.matrix(roundtrip))
})

test_that("sparse matrix saving handles dimnames", {
    X <- Matrix::rsparsematrix(20, 30, 0.2)
    rownames(X) <- LETTERS[1:20]
    colnames(X) <- 1:30

    temp <- saveDelayed(DelayedArray(X))
    roundtrip <- loadDelayed(temp)
    expect_identical(as(X, "SVT_SparseMatrix"), roundtrip@seed)
})

test_that("sparse matrices are saved to external arrays if requested", {
    X <- matrix(rpois(200, 1) - 1L, ncol=20)
    X <- as(X, "SVT_SparseMatrix")
    temp <- saveDelayed(DelayedArray(X), save.external.array=TRUE)
    expect_true(file.exists(file.path(temp, "seeds", 0)))

    roundtrip <- loadDelayed(temp)
    expect_true(is_sparse(roundtrip))
    expect_identical(as.matrix(X), as.matrix(roundtrip))
})

#######################################################
#######################################################

test_that("DelayedAbind works along rows", {
    X <- DelayedArray(matrix(runif(60), ncol=20))
    Y <- DelayedArray(matrix(runif(100), ncol=20))
    Z <- rbind(X, Y)
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    expect_identical(Z, roundtrip)
    expect_s4_class(roundtrip@seed, "DelayedAbind")
})

test_that("DelayedAbind works along columns", {
    A <- DelayedArray(matrix(runif(60), nrow=20))
    B <- DelayedArray(matrix(runif(100), nrow=20))
    C <- DelayedArray(matrix(runif(200), nrow=20)) # throwing in another dataset for some variety.
    Z <- BiocGenerics::cbind(A, B, C)
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    expect_identical(Z, roundtrip)
    expect_s4_class(roundtrip@seed, "DelayedAbind")
})

test_that("DelayedAbind works for 3D arrays", {
    A <- DelayedArray(array(runif(100), c(10, 5, 4)))
    B <- DelayedArray(array(runif(100), c(10, 5, 4)))
    Z <- arbind(A, B)
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    expect_identical(Z, roundtrip)
    expect_s4_class(roundtrip@seed, "DelayedAbind")
})

#######################################################
#######################################################

test_that("DelayedAperm works along rows", {
    X <- DelayedArray(matrix(runif(100), ncol=20))
    Z <- t(X)
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    expect_identical(Z, roundtrip)
    expect_s4_class(roundtrip@seed, "DelayedAperm")
})

test_that("DelayedAperm works for 3D arrays", {
    A <- DelayedArray(array(runif(100), c(10, 5, 4)))
    perm <- c(3L, 1L, 2L)
    Z <- aperm(A, perm)
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    expect_identical(Z, roundtrip)
    expect_s4_class(roundtrip@seed, "DelayedAperm")
})

#######################################################
#######################################################

test_that("DelayedNaryIsoOp works as expected (arithmetic)", {
    X <- DelayedArray(matrix(runif(100), ncol=5))
    Y <- DelayedArray(matrix(runif(100), ncol=5))
    Z <- X * Y
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    expect_identical(Z, roundtrip)
    expect_s4_class(roundtrip@seed, "DelayedNaryIsoOp")
})

test_that("DelayedNaryIsoOp works as expected (comparison)", {
    X <- DelayedArray(matrix(rpois(100, 5), ncol=5))
    Y <- DelayedArray(matrix(rpois(100, 2), ncol=5))
    Z <- X > Y
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    expect_identical(Z, roundtrip)
    expect_s4_class(roundtrip@seed, "DelayedNaryIsoOp")
})

test_that("DelayedNaryIsoOp works as expected (logic)", {
    X <- DelayedArray(matrix(rbinom(100, 1, 0.5) != 0, ncol=5))
    Y <- DelayedArray(matrix(rbinom(100, 1, 0.5) != 0, ncol=5))
    Z <- X | Y
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    expect_identical(Z, roundtrip)
    expect_s4_class(roundtrip@seed, "DelayedNaryIsoOp")
})

test_that("DelayedNaryIsoOp works for 3D arrays", {
    A <- DelayedArray(array(runif(100), c(10, 5, 4)))
    B <- DelayedArray(array(runif(100), c(10, 5, 4)))
    Z <- A <= B
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    expect_identical(Z, roundtrip)
    expect_s4_class(roundtrip@seed, "DelayedNaryIsoOp")
})

test_that("DelayedNaryIsoOp works with multiple seeds", {
    X <- DelayedArray(matrix(rpois(100, 5), ncol=5))
    Y <- DelayedArray(matrix(rpois(100, 2), ncol=5))
    Z <- DelayedArray(matrix(rbinom(100, 1, 0.5) == 0, ncol=5))
    AA <- X - Y + Z
    temp <- saveDelayed(AA)

    roundtrip <- loadDelayed(temp)
    expect_identical(AA, roundtrip)
    expect_s4_class(roundtrip@seed, "DelayedNaryIsoOp")
})

#######################################################
#######################################################

test_that("DelayedSetDimnames works as expected (colnames only)", {
    Z <- DelayedArray(matrix(runif(100), ncol=20))
    colnames(Z) <- LETTERS[1:20]
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    expect_identical(Z, roundtrip)
    expect_s4_class(roundtrip@seed, "DelayedSetDimnames")
})

test_that("DelayedSetDimnames works as expected (rownames only)", {
    Z <- DelayedArray(matrix(runif(100) < 0.1, ncol=20))
    rownames(Z) <- letters[1:5]
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    expect_identical(Z, roundtrip)
    expect_s4_class(roundtrip@seed, "DelayedSetDimnames")
})

test_that("DelayedSetDimnames works as expected (both sets of names)", {
    Z <- DelayedArray(matrix(rpois(100, 5), ncol=20))
    dimnames(Z) <- list(letters[1:5], LETTERS[1:20])
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    expect_identical(Z, roundtrip)
    expect_s4_class(roundtrip@seed, "DelayedSetDimnames")
})

#######################################################
#######################################################

test_that("DelayedSubassign works when all indices are supplied", {
    X <- DelayedArray(matrix(runif(100), ncol=20))
    X[1:2,c(1, 10, 20)] <- matrix(-runif(6), ncol=3)
    temp <- saveDelayed(X)

    roundtrip <- loadDelayed(temp)
    expect_identical(X, roundtrip)
    expect_s4_class(roundtrip@seed, "DelayedSubassign")
})

test_that("DelayedSubassign works when only one index is supplied", {
    X <- DelayedArray(matrix(runif(100), ncol=20))
    X[c(1, 5), ] <- matrix(-rpois(2*ncol(X), 12), ncol=ncol(X))
    temp <- saveDelayed(X)

    roundtrip <- loadDelayed(temp)
    expect_identical(X, roundtrip)
    expect_s4_class(roundtrip@seed, "DelayedSubassign")
})

test_that("DelayedSubassign works when the replacement is a DelayedArray", {
    X <- DelayedArray(matrix(rbinom(100, 1, 0.5) == 0, ncol=20))
    X[1:2,3:5] <- DelayedArray(matrix(-runif(6), ncol=3)) + 1
    temp <- saveDelayed(X)

    roundtrip <- loadDelayed(temp)
    expect_identical(as.matrix(X), as.matrix(roundtrip)) # see comments below for the DelayedUnaryIsoOpStack tests.
    expect_s4_class(roundtrip@seed, "DelayedSubassign")
})

#######################################################
#######################################################

test_that("DelayedSubset works when all indices are supplied", {
    X <- DelayedArray(matrix(runif(100), ncol=20))
    Y <- X[1:2,3:5]
    temp <- saveDelayed(Y)

    roundtrip <- loadDelayed(temp)
    expect_identical(Y, roundtrip)
    expect_s4_class(roundtrip@seed, "DelayedSubset")

    # Trying with out-of-order indices.
    Y <- X[c(1,3,5),c(2,8,6,4)]
    temp <- saveDelayed(Y)
    roundtrip <- loadDelayed(temp)
    expect_identical(Y, roundtrip)
})

test_that("DelayedSubset works when only one index is supplied", {
    X <- DelayedArray(matrix(runif(100), ncol=20))
    Y <- X[,c(20, 1, 3, 15, 9)] 
    temp <- saveDelayed(Y)

    roundtrip <- loadDelayed(temp)
    expect_identical(Y, roundtrip)
    expect_s4_class(roundtrip@seed, "DelayedSubset")
})

#######################################################
#######################################################

# These use as.matrix() in comparisons to avoid problems with environment
# comparisons in the functions stored in OPS.

test_that("DelayedUnaryIsoOpStack works for multiple operations", {
    X <- DelayedArray(matrix(runif(100), ncol=20))
    Z <- abs(X - 0.5)
    expect_s4_class(Z@seed, "DelayedUnaryIsoOpStack")
    expect_type(Z@seed@seed, "double")

    temp <- saveDelayed(Z)
    roundtrip <- loadDelayed(temp)
    expect_identical(as.matrix(Z), as.matrix(roundtrip))
    expect_s4_class(roundtrip@seed, "DelayedUnaryIsoOpStack")
})

test_that("DelayedUnaryIsoOpStack works on both sides", {
    X <- DelayedArray(matrix(rpois(100, 5) + 1L, ncol=20))
    Z <- 5 / X
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    expect_identical(as.matrix(Z), as.matrix(roundtrip))
    expect_s4_class(roundtrip@seed, "DelayedUnaryIsoOpStack")

    # Trying on the other side with another non-commutative op.
    Z <- X - 10
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    expect_identical(as.matrix(Z), as.matrix(roundtrip))
    expect_s4_class(roundtrip@seed, "DelayedUnaryIsoOpStack")
})

test_that("DelayedUnaryIsoOpStack works for log", {
    X <- DelayedArray(matrix(runif(100), ncol=20))
    Z <- log(X)
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    expect_identical(as.matrix(Z), as.matrix(roundtrip))
    expect_s4_class(roundtrip@seed, "DelayedUnaryIsoOpStack")

    # Works with non-default base.
    Z <- log(X, base=3)
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    expect_identical(as.matrix(Z), as.matrix(roundtrip))
    expect_s4_class(roundtrip@seed, "DelayedUnaryIsoOpStack")
})

test_that("DelayedUnaryIsoOpStack works (comparisons)", {
    X <- DelayedArray(matrix(runif(1000), ncol=20))
    Z <- 0.5 < X
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    expect_identical(as.matrix(Z), as.matrix(roundtrip))
    expect_s4_class(roundtrip@seed, "DelayedUnaryIsoOpStack")

    Z <- X <= 0.2
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    expect_identical(as.matrix(Z), as.matrix(roundtrip))
    expect_s4_class(roundtrip@seed, "DelayedUnaryIsoOpStack")
})

test_that("DelayedUnaryIsoOpStack works (logic operations)", {
    X <- DelayedArray(matrix(runif(1000), ncol=20) > 0.5)
    Z <- X & TRUE
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    expect_identical(as.matrix(Z), as.matrix(roundtrip))
    expect_s4_class(roundtrip@seed, "DelayedUnaryIsoOpStack")

    # Same for the ||.
    Z <- X | FALSE
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    expect_identical(as.matrix(Z), as.matrix(roundtrip))
    expect_s4_class(roundtrip@seed, "DelayedUnaryIsoOpStack")

    # Same for !
    Z <- !X
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    expect_identical(as.matrix(Z), as.matrix(roundtrip))
    expect_s4_class(roundtrip@seed, "DelayedUnaryIsoOpStack")
})

test_that("DelayedUnaryIsoOpStack works for unary arithmetic", {
    X <- DelayedArray(matrix(rnorm(1000), ncol=20))
    Z <- -X
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    expect_identical(as.matrix(Z), as.matrix(roundtrip))
    expect_s4_class(roundtrip@seed, "DelayedUnaryIsoOpStack")

    Z <- +X
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    expect_identical(as.matrix(Z), as.matrix(roundtrip))
    expect_type(roundtrip@seed, "double")
})

test_that("DelayedUnaryIsoOpStack works for other unary operations", {
    suppressWarnings(X <- DelayedArray(matrix(log(rnorm(1000)), ncol=20)))
    Z <- is.nan(X)
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    suppressWarnings(expect_identical(as.matrix(Z), as.matrix(roundtrip)))
    expect_s4_class(roundtrip@seed, "DelayedUnaryIsoOpStack")
})

test_that("DelayedUnaryIsoOpStack works for Math2", {
    X <- DelayedArray(matrix(rnorm(1000) * 10, ncol=20))
    Z <- round(X)
    temp <- tempfile(fileext=".h5")
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    expect_identical(as.matrix(Z), as.matrix(roundtrip))
    expect_s4_class(roundtrip@seed, "DelayedUnaryIsoOpStack")

    # Throwing in some non-standard digits.
    Z <- signif(X, digits=3)
    temp <- tempfile(fileext=".h5")
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    expect_identical(as.matrix(Z), as.matrix(roundtrip))
    expect_s4_class(roundtrip@seed, "DelayedUnaryIsoOpStack")
})

#######################################################
#######################################################

test_that("DelayedUnaryIsoOpWithArgs works as expected", {
    X <- DelayedArray(matrix(runif(100), ncol=20))
    Z <- X - runif(5)
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    expect_identical(as.matrix(Z), as.matrix(roundtrip))
    expect_s4_class(roundtrip@seed, "DelayedUnaryIsoOpWithArgs")
})

test_that("DelayedUnaryIsoOpWithArgs works along the other dimension", {
    X <- DelayedArray(matrix(runif(100), ncol=20))
    Z <- X - runif(5)
    temp <- saveDelayed(Z)

    # Manually injecting along=1, because we can't actually seem to stage a
    # DelayedArray directly with along=1; calling DelayedArray::sweep does
    # a double-transpose instead.
    library(rhdf5)
    vec <- runif(20)
    (function() {
        fhandle <- H5Fopen(file.path(temp, "array.h5"))
        on.exit(H5Fclose(fhandle), add=TRUE, after=FALSE)
        ghandle <- H5Gopen(fhandle, "delayed_array")
        on.exit(H5Gclose(ghandle), add=TRUE, after=FALSE)
        H5Ldelete(ghandle, "along")
        h5_write_vector(ghandle, "along", 1, type="H5T_NATIVE_UINT8", scalar=TRUE)
        H5Ldelete(ghandle, "value")
        dhandle <- h5_write_vector(ghandle, "value", vec, type="H5T_NATIVE_DOUBLE", emit=TRUE)
        on.exit(H5Dclose(dhandle), add=TRUE, after=FALSE)
        h5_write_attribute(dhandle, "type", "FLOAT", scalar=TRUE)
    })()

    roundtrip <- loadDelayed(temp)
    expected <- sweep(as.matrix(X), MARGIN=2, STATS=vec, FUN="-")
    expect_identical(as.matrix(expected), as.matrix(roundtrip))
})

test_that("DelayedUnaryIsoOpWithArgs handles logical renaming", {
    X <- DelayedArray(matrix(runif(100) > 0.5, ncol=20))
    Z <- X & runif(5) > 0.5
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    expect_identical(as.matrix(Z), as.matrix(roundtrip))
    expect_s4_class(roundtrip@seed, "DelayedUnaryIsoOpWithArgs")
})

test_that("DelayedUnaryIsoOpWithArgs handles NAs correctly", {
    X <- DelayedArray(matrix(runif(100), ncol=20))
    vec <- runif(5)
    vec[1] <- NA
    Z <- X - vec
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    expect_identical(as.matrix(Z), as.matrix(roundtrip))
    expect_s4_class(roundtrip@seed, "DelayedUnaryIsoOpWithArgs")

    # Works with non-default NAs
    vec <- 1:5
    Z <- X / vec
    temp <- saveDelayed(Z)

    library(rhdf5)
    (function() {
        fhandle <- H5Fopen(file.path(temp, "array.h5"))
        on.exit(H5Fclose(fhandle), add=TRUE)
        dhandle <- H5Dopen(fhandle, "delayed_array/value")
        on.exit(H5Dclose(dhandle), add=TRUE)
        h5_write_attribute(dhandle, "missing_placeholder", 3, type="H5T_NATIVE_DOUBLE", scalar=TRUE)
    })()

    roundtrip <- loadDelayed(temp)
    vec[3] <- NA
    expected <- X / vec
    expect_identical(as.matrix(expected), as.matrix(roundtrip))
    expect_s4_class(roundtrip@seed, "DelayedUnaryIsoOpWithArgs")
})

test_that("DelayedUnaryIsoOpWithArgs handles side-ness correctly", {
    X <- DelayedArray(matrix(runif(100), ncol=20))
    Z <- runif(5) / X
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    expect_identical(as.matrix(Z), as.matrix(roundtrip))
    expect_s4_class(roundtrip@seed, "DelayedUnaryIsoOpWithArgs")
})

test_that("DelayedUnaryIsoOpWithArgs handles repeated operations correctly (same side)", {
    X <- DelayedArray(matrix(runif(100), ncol=20))
    Z <- (X + runif(5)) + runif(5)
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    expect_identical(as.matrix(Z), as.matrix(roundtrip))
    expect_s4_class(roundtrip@seed, "DelayedUnaryIsoOpWithArgs")
})

test_that("DelayedUnaryIsoOpWithArgs handles repeated operations correctly (opposite sides)", {
    X <- DelayedArray(matrix(runif(100), ncol=20))
    Z <- runif(5) + (X + runif(5))
    temp <- saveDelayed(Z)

    roundtrip <- loadDelayed(temp)
    expect_identical(as.matrix(Z), as.matrix(roundtrip))
    expect_s4_class(roundtrip@seed, "DelayedUnaryIsoOpWithArgs")
})

#######################################################
#######################################################

test_that("saving of a LowRankMatrix works correctly", {
    left <- matrix(rnorm(100000), ncol=20)
    right <- matrix(rnorm(50000), ncol=20)

    library(BiocSingular)
    thing <- LowRankMatrix(left, right)
    temp <- saveDelayed(thing)

    out <- loadDelayed(temp)
    expect_identical(thing, out)
})

test_that("saving of a ResidualMatrix works correctly", {
    y <- rsparsematrix(80, 50, 0.5)
    design <- model.matrix(~gl(8, 10))
    thing <- ResidualMatrix::ResidualMatrix(y, design=design)

    # Round-trips properly.
    temp <- saveDelayed(thing)
    out <- loadDelayed(temp)
    out@seed@.matrix <- as(out@seed@.matrix, "dgCMatrix")
    expect_identical(thing, out)
    expect_s4_class(out, "ResidualMatrix")

    # Works when transposed.
    thing2 <- t(thing)
    temp2 <- saveDelayed(thing2)
    out <- loadDelayed(temp2)
    out@seed@.matrix <- as(out@seed@.matrix, "dgCMatrix")
    expect_identical(thing2, out)
    expect_s4_class(out, "ResidualMatrix")

    # Same result if we ignore the type hint.
    (function() {
        fhandle <- H5Fopen(file.path(temp, "array.h5"), "H5F_ACC_RDWR")
        on.exit(H5Fclose(fhandle))
        ghandle <- H5Gopen(fhandle, "delayed_array")
        on.exit(H5Gclose(ghandle))
        H5Ldelete(ghandle, "_r_type_hint")
    })()

    out <- loadDelayed(temp)
    expect_false(is(out, "ResidualMatrix"))
    expect_identical(unname(as.matrix(thing)), unname(as.matrix(out)))

    (function() {
        fhandle <- H5Fopen(file.path(temp2, "array.h5"), "H5F_ACC_RDWR")
        on.exit(H5Fclose(fhandle))
        ghandle <- H5Gopen(fhandle, "delayed_array")
        on.exit(H5Gclose(ghandle))
        H5Ldelete(ghandle, "_r_type_hint")
    })()

    out <- loadDelayed(temp2)
    expect_false(is(out, "ResidualMatrix"))
    expect_identical(unname(as.matrix(thing2)), unname(as.matrix(out)))
})
