#' @importFrom BiocGenerics type
to_array_type <- function(x) {
    switch(type(x),
        integer="integer",
        double="number",
        logical="boolean",
        character="string",
        "other"
    )
}

array_type <- to_array_type

from_array_type <- function(x) {
    switch(x, 
        integer="integer",
        number="double",
        boolean="logical",
        string="character"
    )
}

save_names <- function(handle, x, group = "names", transpose=FALSE) {
    d <- dimnames(x)
    if (is.null(d) || all(vapply(d, is.null, TRUE))) {
        return(NULL)
    }

    if (transpose) { # for the HDF5 array transposition.
        d <- rev(d)
    }

    ghandle <- H5Gcreate(handle, group)
    on.exit(H5Gclose(ghandle), add=TRUE, after=FALSE)
    for (i in seq_along(d)) {
        current <- d[[i]]
        if (!is.null(current)) {
            h5_write_vector(ghandle, as.character(i - 1L), current)
        }
    }
}

load_names <- function(handle, ndim, group = "names") {
    if (!h5_object_exists(handle, group)) {
        return(NULL)
    }

    ghandle <- H5Gopen(handle, group)
    on.exit(H5Gclose(ghandle), add=TRUE, after=FALSE)
    all.named <- h5ls(ghandle, datasetinfo=FALSE, recursive=FALSE)

    names <- vector("list", ndim)
    for (y in all.named$name) {
        names[[as.integer(y) + 1L]] <- h5_read_vector(ghandle, y)
    }

    names
}
