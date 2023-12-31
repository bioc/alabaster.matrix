% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readArray.R
\name{readArray}
\alias{readArray}
\alias{loadArray}
\title{Read a dense array from disk}
\usage{
readArray(path, metadata, array.output.type = NULL, ...)
}
\arguments{
\item{path}{String containing a path to a directory, itself created by the \code{\link{saveObject}} method for a dense array.}

\item{metadata}{Named list of metadata for this object, see \code{\link{readObject}} for more details.}

\item{array.output.type}{String specifying the output type for this function.
This can be \code{"array"} or \code{"ReloadedArray"} (the default).}

\item{...}{Further arguments, ignored.}
}
\value{
A multi-dimensional array-like object.
}
\description{
Read a dense high-dimensional array from its on-disk representation.
This is usually not directly called by users, but is instead called by dispatch in \code{\link{readObject}}.
}
\details{
By default, a file-backed \linkS4class{ReloadedArray} is returned to save memory and to preserve the provenance of the object.
Users can set \code{array.output.type="array"} to force \code{readArray} to load all data into memory and return an ordinary R array instead.
}
\examples{
arr <- array(rpois(10000, 10), c(50, 20, 10))
dimnames(arr) <- list(
   paste0("GENE_", seq_len(nrow(arr))),
   letters[1:20],
   NULL
)

dir <- tempfile()
saveObject(arr, dir)
readObject(dir)

}
\seealso{
\code{"\link{saveObject,array-method}"}, to create the directory and its contents.
}
\author{
Aaron Lun
}
