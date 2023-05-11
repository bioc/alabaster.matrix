#' @export
#' @import methods
setClass("WrapperArraySeed", contains="VIRTUAL", slots=c(seed="ANY"))

#' @export
#' @importClassesFrom DelayedArray DelayedAbind
setClass("AmalgamatedArraySeed", contains="DelayedAbind", slots=c(samples = "character"))

#' @export
#' @importClassesFrom DelayedArray DelayedArray
setClass("AmalgamatedArray", contains="DelayedArray", slots=c(seed = "AmalgamatedArraySeed"))

#' @export
#' @importClassesFrom DelayedArray DelayedMatrix
setClass("AmalgamatedMatrix", contains=c("AmalgamatedArray", "DelayedMatrix"))
