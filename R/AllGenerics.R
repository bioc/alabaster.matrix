#' @export
setGeneric("loadWrapperArray", function(meta, project) standardGeneric("loadWrapperArray"), signature="project")

#' @export
#' @rdname storeDelayedObject
setGeneric("storeDelayedObject", function(x, handle, name, ...) standardGeneric("storeDelayedObject"))
