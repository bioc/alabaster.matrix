# Sticking this test here because of how the dependencies work.
# library(testthat); library(alabaster.matrix); source("test-stage-df.R")

library(alabaster.base)
library(S4Vectors)

test_that("staging of arrays within DFs works correctly", {
    tmp <- tempfile()
    dir.create(tmp)

    input <- DataFrame(A=1:3, X=0, Y=0)
    input$X <- array(runif(3)) # 1D array
    input$Y <- cbind(runif(3)) # matrix with 1 column 
    input$Z <- cbind(V=runif(3), Z=rnorm(3)) # matrix with 2 columns.
    info <- stageObject(input, tmp, path="WHEE")

    roundtrip <- loadDataFrame(info, project=tmp)
    input$X <- as.numeric(input$X) # strip 1D arrays as these are rarely desirable.
    roundtrip$Y <- as.matrix(roundtrip$Y)
    roundtrip$Z <- as.matrix(roundtrip$Z)
    expect_equal(roundtrip, input)
})

test_that("staging of arrays continues to work with character matrices", {
    tmp <- tempfile()
    dir.create(tmp)

    input <- DataFrame(A=1:3, X=0, Y=0)
    input$X <- array(letters[1:3]) # 1D array
    input$Y <- cbind(LETTERS[1:3]) # matrix with 1 column 
    input$Z <- cbind(V=letters[4:6], Z=LETTERS[4:6]) # matrix with 2 columns.
    info <- stageObject(input, tmp, path="WHEE")

    roundtrip <- loadDataFrame(info, project=tmp)
    input$X <- as.character(input$X) # strip 1D arrays as these are rarely desirable.
    roundtrip$Y <- as.matrix(roundtrip$Y)
    roundtrip$Z <- as.matrix(roundtrip$Z)
    expect_equal(roundtrip, input)
})

