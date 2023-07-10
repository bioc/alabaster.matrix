#' Amalgamated array class
#'
#' Implements an amalgamated array, equivalent to a delayed combination of DelayedArray objects.
#' It allows \code{\link{stageObject}} to save a combination of multiple matrices without actually aggregating their data into a single file.
#'
#' @section Constructors:
#' \code{AmalgamatedArraySeed(..., along=1)} accepts any number of named array-like objects and returns a AmalgamatedArraySeed.
#' Each object corresponds to a block and should be named accordingly; names should be unique and non-empty.
#' The \code{along} argument specifies the dimension in which matrices should be combined - for matrices, this is 1 for rows, 2 for columns.
#'
#' \code{AmalgamatedArray(..., along=1)} accepts any number of named array-like objects and returns a AmalgamatedArray.
#' Alternatively, a single AmalgamatedArraySeed may be provided in \code{...}. 
#' 
#' @section Functions:
#' \code{componentNames(x)} will return a character vector of names of component arrays in a AmalgamatedArray(Seed) object \code{x}.
#'
#' \code{extractComponents(x)} will return a named list of array-like objects,
#' corresponding to the component arrays used to construct the AmalgamatedArray(Seed) object \code{x}.
#'
#' \code{\link{stageObject}(x, dir, path, child = FALSE)} will save the AmalgamatedArray \code{x} and its components into the \code{path} inside \code{dir}.
#' Each component array is staged into its own subdirectory inside \code{path}.
#'
#' @section Comments on usage:
#' The AmalgamatedArraySeed is closely related to (and in fact, is a subclass of) the \linkS4class{DelayedAbind} class.
#' This means that we can leverage many of the \pkg{DelayedArray} methods for handling the delayed bind.
#' In theory, we could just use a DelayedAbind directly and save it with \pkg{chihaya} in \code{\link{stageObject}} (via \code{\link{preserveDelayedOperations}(TRUE)}).
#' However, this provides fewer opportunities for tracking and manipulating the samples.
#' It also saves the per-sample matrices into a single file, which eliminates possibilities for per-file deduplication and linking, e.g., with \code{\link{recycleHdf5Files}(TRUE)}.
#' 
#' @author Aaron Lun
#' @examples
#' first <- Matrix::rsparsematrix(10, 10, 0.1)
#' second <- Matrix::rsparsematrix(10, 20, 0.1)
#' mat <- AmalgamatedArray(list(foo = first, bar = second), along=2)
#' mat
#'
#' componentNames(mat)
#' out <- extractComponents(mat)
#' lapply(out, dim)
#'
#' @aliases
#' AmalgamatedArray
#' AmalgamatedArray-class
#' AmalgamatedMatrix-class
#' AmalgamatedArraySeed
#' AmalgamatedArraySeed-class
#' DelayedArray,AmalgamatedArraySeed-method
#' matrixClass,AmalgamatedArray-method
#' componentNames
#' extractComponents
#' stageObject,AmalgamatedArray-method
#' coerce,AmalgamatedArray,AmalgamatedMatrix-method
#' coerce,AmalgamatedMatrix,AmalgamatedArray-method
#' 
#' @name AmalgamatedArray
NULL

#' @export
#' @importFrom DelayedArray arbind acbind
AmalgamatedArraySeed <- function(components, along = 1) {
    sample.names <- names(components)
    if (anyDuplicated(sample.names) || any(sample.names == "")) {
        stop("sample names must be unique and non-empty in a AmalgamatedArraySeed")
    }

    FUN <- if (along == 1) arbind else acbind
    combined <- do.call(FUN, lapply(components, DelayedArray))
    new("AmalgamatedArraySeed", combined@seed, samples = sample.names)
}

#' @export
#' @importFrom DelayedArray DelayedArray new_DelayedArray
setMethod("DelayedArray", "AmalgamatedArraySeed", function(seed) new_DelayedArray(seed, Class="AmalgamatedArray"))

#' @export
#' @importFrom DelayedArray matrixClass
setMethod("matrixClass", "AmalgamatedArray", function(x) "AmalgamatedMatrix")

# Overrides copied from DelayedArray::ConstantArray.
#' @importFrom S4Vectors new2
setAs("AmalgamatedArray", "AmalgamatedMatrix", function(from) new2("AmalgamatedMatrix", from))
setAs("AmalgamatedMatrix", "AmalgamatedArray", function(from) from)

#' @export
AmalgamatedArray <- function(components, along = 1) {
    if (!is(components, "AmalgamatedArraySeed")) {
        seed <- AmalgamatedArraySeed(components, along = along)
    } else {
        seed <- components
    }
    DelayedArray(seed)
}

#' @export
componentNames <- function(x) {
    if (is(x, "DelayedMatrix")) {
        x <- x@seed
    }
    x@samples
}

#' @export
extractComponents <- function(x) {
    if (is(x, "DelayedMatrix")) {
        x <- x@seed
    }
    out <- x@seeds
    names(out) <- x@samples
    out
}

#' @export
#' @importFrom alabaster.base .stageObject .writeMetadata
#' @importFrom DelayedArray DelayedArray
setMethod("stageObject", "AmalgamatedArray", function(x, dir, path, child = FALSE) {
    dir.create(file.path(dir, path), showWarnings=FALSE)

    seed <- x@seed
    seeds <- seed@seeds
    components <- vector("list", length(seeds))

    for (i in seq_along(seeds)) {
        current.seed <- seeds[[i]]
        out <- try(current.seed[0,0], silent=TRUE)
        if (is(out, "try-error")) {
            current.seed <- DelayedArray(current.seed)
        }
        meta <- .stageObject(current.seed, dir, paste0(path, "/component", i), child = TRUE)
        components[[i]] <- list(name = seed@samples[i], resource = .writeMetadata(meta, dir))
    }

    list(
        `$schema` = "amalgamated_array/v1.json",
        path = paste0(path, "/array.json"),
        is_child = child,
        `array` = .grab_array_type(x),
        amalgamated_array = list(
            along = x@seed@along - 1L,
            components = components
        )
    )
})

#' @importFrom alabaster.base acquireMetadata .loadObject
.load_amalgamated_array <- function(info, project) {
    comp <- info$amalgamated_array$components
    all.names <- character(length(comp))

    for (i in seq_along(comp)) {
        current <- comp[[i]]
        imeta <- acquireMetadata(project, current$resource$path)
        all.names[i] <- current$name
        comp[[i]] <- .loadObject(imeta, project)
    }

    names(comp) <- all.names
    AmalgamatedArraySeed(comp, along = info$amalgamated_array$along + 1L)
}

