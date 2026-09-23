#' @export
#' @importClassesFrom DelayedArray DelayedAbind
setClass("AmalgamatedArraySeed", contains="DelayedAbind", slots=c(samples = "character"))

#' @export
#' @importClassesFrom DelayedArray DelayedArray
setClass("AmalgamatedArray", contains="DelayedArray", slots=c(seed = "AmalgamatedArraySeed"))

#' @export
#' @importClassesFrom DelayedArray DelayedMatrix
setClass("AmalgamatedMatrix", contains=c("AmalgamatedArray", "DelayedMatrix"))

#' Delayed no-op base class 
#'
#' A DelayedArray layer that does not modify the contents of its seed.
#' Generics like \code{\link[DelayedArray]{extract_array}} simply call the corresponding method for the seed.
#' DelayedNoOp subclasses are used by \pkg{alabaster.matrix} to provide metadata that persists throughout delayed operations,
#' see \linkS4class{ReloadedArraySeed} and \linkS4class{ToBeRealizedArraySeed} for more details.
#'
#' We expose this subclass to facilitate traversal of a DelayedArray's operation tree that contains \pkg{alabaster.matrix} classes.
#' Specifically, these no-op instances can be easily detected and skipped without having to list each of the individual subclasses.
#' For example, packages like \pkg{beachmat} can check if the current delayed layer is a no-op class and immediately move onto processing the seed.
#'
#' @export
#' @importClassesFrom DelayedArray DelayedUnaryIsoOp
#' @aliases DelayedNoOp-class
#' @name DelayedNoOp
setClass("DelayedNoOp", contains="DelayedUnaryIsoOp")

#' @export
setClass("ReloadedArraySeed", contains="DelayedNoOp", slots=c(path="character"))

#' @export
setClass("ReloadedArray", contains="DelayedArray", slots=c(seed="ReloadedArraySeed"))

#' @export
#' @importClassesFrom DelayedArray DelayedMatrix
setClass("ReloadedMatrix", contains=c("ReloadedArray", "DelayedMatrix"))

#' @export
setClass("ToBeRealizedArraySeed", contains="DelayedNoOp")

#' @export
setClass("ToBeRealizedArray", contains="DelayedArray", slots=c(seed="ToBeRealizedArraySeed"))

#' @export
#' @importClassesFrom DelayedArray DelayedMatrix
setClass("ToBeRealizedMatrix", contains=c("ToBeRealizedArray", "DelayedMatrix"))
