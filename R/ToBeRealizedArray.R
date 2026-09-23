#' To-be-realized array 
#'
#' A DelayedArray layer that indicates that its seed is to be realized by \code{\link{saveObject,DelayedArray-method}}.
#' This is intended for fine-grained control over which delayed operations are to be preserved in an object with multiple arrays, e.g., a SummarizedExperiment.
#'
#' @param seed Seed representing an array-like object, typically a \link[DelayedArray]{DelayedArray}.
#' @param x A \linkS4class{ToBeRealizedArray} object.
#' @param path String containing a path to a directory in which to save \code{x}.
#' @param array.dedup.session,array.dedup.action Arguments controlling deduplication of \code{x}, see \code{?"\link{saveObject,array-method}"} for details.
#' @param ... Further arguments to pass to the relevant \code{\link{saveObject}} methods for saving dense arrays and sparse matrices. 
#'
#' @return
#' For the constructors, an instance of the \linkS4class{ToBeRealizedArraySeed} or \linkS4class{ToBeRealizedArray}.
#'
#' For \code{saveObject}, \code{x} is saved to \code{path} and \code{NULL} is returned. 
#'
#' @details
#' The ToBeRealizedArraySeed is a \linkS4class{DelayedNoOp} subclass that will just forward all operations to the underlying \code{seed}.
#' Its purpose is to indicate to \code{\link{saveObject}} that any delayed operations in its \code{seed} should not be preserved,
#' even if \code{DelayedArray.preserve.ops=TRUE} in \code{\link{saveObject,DelayedArray-method}}.
#' This is occasionally useful when saving a complicated object with many arrays like a SummarizedExperiment,
#' where we want to preserve delayed operations in some arrays while realizing operations in others.
#'
#' @author Aaron Lun
#'
#' @examples
#' library(DelayedArray)
#' x1 <- matrix(runif(100), 5, 20)
#' x2 <- matrix(runif(250), 5, 50)
#' x <- cbind(DelayedArray(x1), DelayedArray(x2))
#' y0 <- matrix(runif(300), 10, 30)
#' y <- log1p(x)
#' big.object <- list(foo=x, bar=y)
#'
#' # Imagine we want to save 'big.object', preserving the delayed log1p
#' # for 'y' but realizing the delayed cbind for 'x'. This would not be
#' # possible with just DelayedArray.preserve.ops=TRUE as both sets of
#' # delayed operations would be preserved:
#' tmp <- tempfile()
#' saveObject(big.object, tmp, DelayedArray.preserve.ops=TRUE)
#' list.files(tmp, recursive=TRUE)
#' readObjectFile(file.path(tmp, "other_contents", "0"))$type
#' readObjectFile(file.path(tmp, "other_contents", "1"))$type
#'
#' # Instead we wrap 'x' in a ToBeRealizedArray, which forces its realization.
#' big.object2 <- list(foo=ToBeRealizedArray(x), bar=y)
#' tmp2 <- tempfile()
#' saveObject(big.object2, tmp2, DelayedArray.preserve.ops=TRUE)
#' list.files(tmp2, recursive=TRUE)
#' readObjectFile(file.path(tmp2, "other_contents", "0"))$type
#' readObjectFile(file.path(tmp2, "other_contents", "1"))$type
#'
#' @aliases ToBeRealizedArraySeed-class
#' @aliases ToBeRealizedArray-class
#' @aliases ToBeRealizedMatrix-class
#' @aliases DelayedArray,ToBeRealizedArraySeed-method
#' @aliases matrixClass,ToBeRealizedArray-method
#' @aliases coerce,ToBeRealizedArray,ToBeRealizedMatrix-method
#' @aliases coerce,ToBeRealizedMatrix,ToBeRealizedArray-method
#' @aliases saveObject,ToBeRealizedArray-method
#'
#' @export
ToBeRealizedArraySeed <- function(seed) {
    if (is(seed, "ToBeRealizedArraySeed")) {
        return(seed)
    }
    while (is(seed, "DelayedArray")) {
        seed <- seed@seed
    }
    new("ToBeRealizedArraySeed", seed=seed)
}

#' @export
#' @rdname ToBeRealizedArraySeed
ToBeRealizedArray <- function(seed) {
    DelayedArray(ToBeRealizedArraySeed(seed))
}

#' @export
#' @importFrom DelayedArray DelayedArray new_DelayedArray
setMethod("DelayedArray", "ToBeRealizedArraySeed", function(seed) new_DelayedArray(seed, Class="ToBeRealizedArray"))

#' @export
#' @importFrom DelayedArray matrixClass
setMethod("matrixClass", "ToBeRealizedArray", function(x) "ToBeRealizedMatrix")

# Overrides copied from DelayedArray::ConstantArray.
#' @importFrom S4Vectors new2
setAs("ToBeRealizedArray", "ToBeRealizedMatrix", function(from) new2("ToBeRealizedMatrix", from))
setAs("ToBeRealizedMatrix", "ToBeRealizedArray", function(from) from)

#' @export
#' @rdname ToBeRealizedArraySeed
setMethod("saveObject", "ToBeRealizedArray", function(x, path, array.dedup.session=NULL, array.dedup.action=NULL, ...) {
    if (is_sparse(x)) {
        .save_compressed_sparse_matrix(x, path, array.dedup.session=array.dedup.session, array.dedup.action=array.dedup.action, ...)
    } else {
        .save_array(x, path, array.dedup.session=array.dedup.session, array.dedup.action=array.dedup.action, ...)
    }
    invisible(NULL)
})
