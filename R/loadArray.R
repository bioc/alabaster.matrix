#' Load high-dimensional arrays
#'
#' Default loading of arrays from on-disk formats, using the corresponding \code{\link{stageObject}} method. 
#' It should not be necessary for users to call this function manually. 
#'
#' @param info Named list containing the metadata for this array.
#' @param project Any argument accepted by the acquisition functions, see \code{?\link{acquireFile}}. 
#' By default, this should be a string containing the path to a staging directory.
#' 
#' @return A multi-dimensional object (usually a \linkS4class{DelayedMatrix}) containing the array data.
#'
#' @author Aaron Lun
#'
#' @examples
#' dir <- tempfile()
#' dir.create(dir)
#'
#' arr <- array(rpois(10000, 10), c(50, 20, 10))
#' dimnames(arr) <- list(
#'    paste0("GENE_", seq_len(nrow(arr))),
#'    letters[1:20],
#'    NULL
#' )
#'
#' path <- "whee"
#' info <- stageObject(arr, dir, path)
#' loadArray(info, project=dir)
#' 
#' @export
#' @importFrom DelayedArray DelayedArray
#' @importFrom alabaster.base acquireFile
loadArray <- function(info, project) {
    path <- acquireFile(project, info$path)
    seed <- .createRawArraySeed(info, path, names=TRUE)
    DelayedArray(seed)
}
