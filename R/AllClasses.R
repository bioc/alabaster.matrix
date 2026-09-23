#' @export
#' @importClassesFrom DelayedArray DelayedAbind
setClass("AmalgamatedArraySeed", contains="DelayedAbind", slots=c(samples = "character"))

#' @export
#' @importClassesFrom DelayedArray DelayedArray
setClass("AmalgamatedArray", contains="DelayedArray", slots=c(seed = "AmalgamatedArraySeed"))

#' @export
#' @importClassesFrom DelayedArray DelayedMatrix
setClass("AmalgamatedMatrix", contains=c("AmalgamatedArray", "DelayedMatrix"))

#' @export
setClass("ReloadedArraySeed", contains="DelayedUnaryIsoOp", slots=c(path="character"))

#' @export
setClass("ReloadedArray", contains="DelayedArray", slots=c(seed="ReloadedArraySeed"))

#' @export
#' @importClassesFrom DelayedArray DelayedMatrix
setClass("ReloadedMatrix", contains=c("ReloadedArray", "DelayedMatrix"))

#' @export
setClass("ToBeRealizedArraySeed", contains="DelayedUnaryIsoOp", slots=c(path="character"))

#' @export
setClass("ToBeRealizedArray", contains="DelayedArray", slots=c(seed="ToBeRealizedArraySeed"))

#' @export
#' @importClassesFrom DelayedArray DelayedMatrix
setClass("ToBeRealizedMatrix", contains=c("ToBeRealizedArray", "DelayedMatrix"))
