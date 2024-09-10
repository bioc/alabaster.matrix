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

path_components <- function(path) {
    output <- character()
    while (TRUE) {
        base <- basename(path)
        dpath <- dirname(path)
        output <- c(base, output)
        if (dpath == path) {
            break
        }
        path <- dpath
    }
    output
}

clone_duplicate <- function(src, dest, action) {
    dir.create(dest)
    manifest <- list.files(src, recursive=TRUE)

    if (action == "relsymlink") {
        # Find relative path from one to the other.
        src <- normalizePath(src, mustWork=TRUE)
        src.comp <- path_components(src)
        src.len <- length(src.comp)

        dest <- normalizePath(dest, mustWork=TRUE)
        dest.comp <- path_components(dest)
        dest.len <- length(dest.comp)

        counter <- 0L
        for (i in seq_len(min(src.len, dest.len))) {
            if (src.comp[i] != dest.comp[i]) {
                counter <- i - 1L
                break
            }
        }

        base2base <- do.call(file.path, as.list(c(rep("..", dest.len - counter), src.comp[(counter+1):src.len])))
        pwd <- getwd()
        on.exit(setwd(pwd), add=TRUE, after=FALSE)
        setwd(dest)

        for (y in manifest) {
            if (!file.symlink(file.path(base2base, y), y)) {
                stop("failed to link '", y, "' from '", src, "' to '", dest, "'")
            }
        }
        return(NULL)
    }

    if (action == "link") {
        fun <- function(from, to) file.link(from, to) || file.copy(from, to)
        msg <- "copy or link"
    } else if (action == "copy") {
        fun <- file.copy
        msg <- "copy"
    } else if (action == "symlink") {
        fun <- file.symlink
        msg <- "link"
    }

    for (y in manifest) {
        if (!fun(file.path(src, y), file.path(dest, y))) {
            stop("failed to ", msg, " '", y, "' from '", src, "' to '", dest, "'")
        }
    }
}

#' @importFrom alabaster.base saveObject
try_altSaveObject <- function(x, ...) {
    if (is.null(altSaveObjectFunction())) {
        fun <- selectMethod(saveObject, signature=class(x), optional=TRUE)
        if (!is.null(fun)) {
            fun(x, ...)
            return(TRUE)
        }
        return(FALSE)
    }

    # We can't use selectMethod() here as we don't know what generic system
    # is being used by altSaveObject, so we just have to try it out.
    fail <- try(altSaveObject(x, ...), silent=TRUE)
    if (!is(fail, "try-error")) {
        return(TRUE)
    }

    # If the failure wasn't due to S4 dispatch, we emit a warning.
    msg <- attr(fail, "condition")
    if (!startsWith(msg$message, "unable to find an inherited method")) {
        warning(msg$message)
    }
    return(FALSE)
}
