optimize_storage <- function(x) {
    tt <- type(x)
    if (tt == "character") {
        optimize_string_storage(x)
    } else if (tt == "double") {
        optimize_float_storage(x)
    } else if (tt == "integer") {
        optimize_integer_storage(x)
    } else if (tt == "logical") {
        optimize_boolean_storage(x)
    } else {
        stop("unsupported type '", tt, "'")
    }
}

###################################################
###################################################

aggregate_range <- function(collated, name) {
    c(
        min(unlist(lapply(collated, function(y) y[[name]][1]))),
        max(unlist(lapply(collated, function(y) y[[name]][2])))
    )
}

aggregate_any <- function(collated, name) {
    any(vapply(collated, function(y) y[[name]], TRUE))
}

aggregate_max <- function(collated, name) {
    max(unlist(lapply(collated, function(y) y[[name]])), na.rm=TRUE)
}

aggregate_sum <- function(collated, name) {
    sum(vapply(collated, function(y) y[[name]], 0L))
}

collect_from_SVT <- function(svt, fun, tt) {
    if (is.null(svt)) {
        attrs <- fun(tt(0))
        attrs$non_zero <- 0L
        return(list(attrs))
    }

    output <- vector("list", length(svt))
    for (i in seq_along(svt)) {
        node <- svt[[i]] 
        if (is.null(node)) {
            val <- tt(0)
        } else {
            val <- node[[2]]
        }
        attrs <- fun(val)
        attrs$non_zero <- length(val)
        output[[i]] <- attrs
    }
    output
}

###################################################
###################################################

setGeneric("collect_integer_attributes", function(x) standardGeneric("collect_integer_attributes"))

.simple_integer_collector <- function(x) {
    list(
        range=suppressWarnings(range(x, na.rm=TRUE)),
        missing=anyNA(x)
    )
}

.combine_integer_attributes <- function(collated) {
    list(
        range=aggregate_range(collated, "range"),
        missing=aggregate_any(collated, "missing")
    )
}

setMethod("collect_integer_attributes", "array", .simple_integer_collector)

setMethod("collect_integer_attributes", "SVT_SparseMatrix", function(x) {
    collated <- collect_from_SVT(x@SVT, .simple_integer_collector, integer)
    output <- list(non_zero=aggregate_sum(collated, "non_zero"))
    c(output, .combine_integer_attributes(collated))
})

#' @importFrom S4Arrays is_sparse
#' @importFrom SparseArray nzdata
setMethod("collect_integer_attributes", "ANY", function(x) {
    output <- list()

    if (is_sparse(x)) {
        collated <- blockApply(x, function(y) {
            out <- .simple_integer_collector(nzdata(y))
            out$non_zero <- length(nzdata(y))
            out
        }, as.sparse=TRUE)
        output$non_zero <- aggregate_sum(collated, "non_zero")
    } else {
        collated <- blockApply(x, .simple_integer_collector)
    }

    c(output, .combine_integer_attributes(collated))
})

optimize_integer_storage <- function(x) {
    attr <- collect_integer_attributes(x)

    if (attr$missing) {
        lower <- attr$range[1]
        upper <- attr$range[2]
        if (is.infinite(lower)) {
            return(list(type="H5T_NATIVE_INT8", placeholder=as.integer(-2^7), size=attr$non_zero))
        }

        if (lower < 0L) {
            if (lower > -2^7 && upper < 2^7) {
                return(list(type="H5T_NATIVE_INT8", placeholder=as.integer(-2^7), size=attr$non_zero))
            } else if (lower > -2^15 && upper < 2^15) {
                return(list(type="H5T_NATIVE_INT16", placeholder=as.integer(-2^15), size=attr$non_zero))
            }
        } else {
            if (upper < 2^8 - 1) {
                return(list(type="H5T_NATIVE_UINT8", placeholder=as.integer(2^8-1), size=attr$non_zero))
            } else if (upper < 2^16 - 1) {
                return(list(type="H5T_NATIVE_UINT16", placeholder=as.integer(2^16-1), size=attr$non_zero))
            }
        }

        return(list(type="H5T_NATIVE_INT32", placeholder=NA_integer_, size=attr$non_zero))

    } else {
        lower <- attr$range[1]
        upper <- attr$range[2]
        if (is.infinite(lower)) {
            return(list(type="H5T_NATIVE_INT8", placeholder=NULL, size=attr$non_zero))
        }

        if (lower < 0L) {
            if (lower >= -2^7 && upper < 2^7) {
                return(list(type="H5T_NATIVE_INT8", placeholder=NULL, size=attr$non_zero))
            } else if (lower >= -2^15 && upper < 2^15) {
                return(list(type="H5T_NATIVE_INT16", placeholder=NULL, size=attr$non_zero))
            }
        } else {
            if (upper < 2^8) {
                return(list(type="H5T_NATIVE_UINT8", placeholder=NULL, size=attr$non_zero))
            } else if (upper < 2^16) {
                return(list(type="H5T_NATIVE_UINT16", placeholder=NULL, size=attr$non_zero))
            }
        }

        return(list(type="H5T_NATIVE_INT32", placeholder=NULL, size=attr$non_zero))
    }
}

###################################################
###################################################

setGeneric("collect_float_attributes", function(x) standardGeneric("collect_float_attributes"))

setMethod("collect_float_attributes", "array", collect_double_attributes)

setMethod("collect_float_attributes", "ddenseMatrix", function(x) collect_double_attributes(x@x))

setMethod("collect_float_attributes", "dsparseMatrix", function(x) {
    out <- collect_double_attributes(x@x)
    out$non_zero <- length(x@x)
    out
})

.combine_float_attributes <- function(collated) {
    output <- list(range=aggregate_range(collated, "range"))
    for (n in c("missing", "non_integer", "has_NaN", "has_Inf", "has_nInf", "has_lowest", "has_highest")) {
        output[[n]] <- aggregate_any(collated, n)
    }
    output
}

setMethod("collect_float_attributes", "SVT_SparseMatrix", function(x) {
    collated <- collect_from_SVT(x@SVT, collect_double_attributes, double)
    output <- list(non_zero=aggregate_sum(collated, "non_zero"))
    c(output, .combine_float_attributes(collated))
})

#' @importFrom S4Arrays is_sparse
#' @importFrom SparseArray nzdata
setMethod("collect_float_attributes", "ANY", function(x) {
    output <- list()
    if (is_sparse(x)) {
        collated <- blockApply(x, function(y) {
            nzd <- nzdata(y)
            out <- collect_double_attributes(nzd)
            out$non_zero <- length(nzd)
            out
        }, as.sparse=TRUE)
        output$non_zero <- aggregate_sum(collated, "non_zero")
    } else {
        collated <- blockApply(x, collect_double_attributes)
    }
    c(output, .combine_float_attributes(collated))
})

optimize_float_storage <- function(x) {
    attr <- collect_float_attributes(x)

    if (attr$missing) {
        if (!attr$non_integer) {
            lower <- attr$range[1]
            upper <- attr$range[2]
            if (lower < 0L) {
                if (lower > -2^7 && upper < 2^7) {
                    return(list(type="H5T_NATIVE_INT8", placeholder=-2^7, size=attr$non_zero))
                } else if (lower > -2^15 && upper < 2^15) {
                    return(list(type="H5T_NATIVE_INT16", placeholder=-2^15, size=attr$non_zero))
                } else if (lower > -2^31 && upper < 2^31) {
                    return(list(type="H5T_NATIVE_INT32", placeholder=-2^31, size=attr$non_zero))
                }
            } else {
                if (upper < 2^8-1) {
                    return(list(type="H5T_NATIVE_UINT8", placeholder=2^8-1, size=attr$non_zero))
                } else if (upper < 2^16-1) {
                    return(list(type="H5T_NATIVE_UINT16", placeholder=2^16-1, size=attr$non_zero))
                } else if (upper < 2^32-1) {
                    return(list(type="H5T_NATIVE_UINT32", placeholder=2^32-1, size=attr$non_zero))
                }
            }
        }

        placeholder <- NULL
        if (!attr$has_NaN) {
            placeholder <- NaN
        } else if (!attr$has_Inf) {
            placeholder <- Inf
        } else if (!attr$has_nInf) {
            placeholder <- -Inf
        } else if (!attr$has_lowest) {
            placeholder <- lowest_double()
        } else if (!attr$has_highest) {
            placeholder <- highest_double()
        }

        # Fallback that just goes through and pulls out all unique values.
        if (is.null(placeholder)) {
            if (is_sparse(x)) {
                u <- Reduce(union, blockApply(x, function(y) unique(nzdata(y))))
            } else {
                u <- Reduce(union, blockApply(x, function(y) unique(as.vector(y))))
            }
            placeholder <- chooseMissingPlaceholderForHdf5(u)
        }

        return(list(type="H5T_NATIVE_DOUBLE", placeholder=placeholder, size=attr$non_zero))

    } else {
        if (!attr$non_integer) {
            lower <- attr$range[1]
            upper <- attr$range[2]
            if (lower < 0L) {
                if (lower >= -2^7 && upper < 2^7) {
                    return(list(type="H5T_NATIVE_INT8", placeholder=NULL, size=attr$non_zero))
                } else if (lower >= -2^15 && upper < 2^15) {
                    return(list(type="H5T_NATIVE_INT16", placeholder=NULL, size=attr$non_zero))
                } else if (lower >= -2^31 && upper < 2^31) {
                    return(list(type="H5T_NATIVE_INT32", placeholder=NULL, size=attr$non_zero))
                }
            } else {
                if (upper < 2^8) {
                    return(list(type="H5T_NATIVE_UINT8", placeholder=NULL, size=attr$non_zero))
                } else if (upper < 2^16) {
                    return(list(type="H5T_NATIVE_UINT16", placeholder=NULL, size=attr$non_zero))
                } else if (upper < 2^32) {
                    return(list(type="H5T_NATIVE_UINT32", placeholder=NULL, size=attr$non_zero))
                }
            }
        }

        return(list(type="H5T_NATIVE_DOUBLE", placeholder=NULL, size=attr$non_zero))
    }
}

###################################################
###################################################

setGeneric("collect_string_attributes", function(x) standardGeneric("collect_string_attributes"))

setMethod("collect_string_attributes", "ANY", function(x) {
    collected <- blockApply(x, function(y) {
        list(
            has_na1=any(y == "NA", na.rm=TRUE),
            has_na2=any(y == "_NA", na.rm=TRUE),
            max_len=suppressWarnings(max(nchar(y, "bytes"), na.rm=TRUE)),
            missing=anyNA(y),
            encoding=unique(Encoding(y))
        )
    })

    list(
        has_na1=aggregate_any(collected, "has_na1"),
        has_na2=aggregate_any(collected, "has_na2"),
        max_len=aggregate_max(collected, "max_len"),
        missing=aggregate_any(collected, "missing"),
        encoding=Reduce(union, lapply(collected, function(y) y$encoding))
    )
})

optimize_string_storage <- function(x) {
    attr <- collect_string_attributes(x)

    placeholder <- NULL
    if (attr$missing) {
        if (!attr$has_na1) {
            placeholder <- "NA"
        } else if (!attr$has_na2) {
            placeholder <- "_NA"
        } else {
            u <- Reduce(union, blockApply(x, function(y) unique(as.vector(y))))
            placeholder <- chooseMissingPlaceholderForHdf5(u)
        }
        attr$max_len <- max(attr$max_len, nchar(placeholder, "bytes"))
    }

    tid <- H5Tcopy("H5T_C_S1")
    H5Tset_strpad(tid, strpad = "NULLPAD")
    H5Tset_size(tid, max(1L, attr$max_len))
    if ("UTF-8" %in% attr$encoding) {
        H5Tset_cset(tid, "UTF8")
    } else {
        H5Tset_cset(tid, "ASCII")
    }

    list(type=tid, placeholder=placeholder)
}

###################################################
###################################################

setGeneric("collect_boolean_attributes", function(x) standardGeneric("collect_boolean_attributes"))

setMethod("collect_boolean_attributes", "ANY", function(x) {
    output <- list()
    if (is_sparse(x)) {
        collated <- blockApply(x, function(x) list(missing=anyNA(nzdata(x)), non_zero=length(nzdata(x))), as.sparse=TRUE)
        output$non_zero <- aggregate_sum(collated, "non_zero")
    } else {
        collated <- list(list(missing=anyNA(x)))
    }
    output$missing <- aggregate_any(collated, "missing")
    output
})

setMethod("collect_boolean_attributes", "lsparseMatrix", function(x) {
    list(missing=anyNA(x), non_zero=length(x@x))
})

setMethod("collect_boolean_attributes", "SVT_SparseMatrix", function(x) {
    collated <- collect_from_SVT(x@SVT, function(vals) { list(missing=anyNA(vals)) }, logical)
    list(
        missing=aggregate_any(collated, "missing"),
        non_zero=aggregate_sum(collated, "non_zero")
    )
})

optimize_boolean_storage <- function(x) {
    attr <- collect_boolean_attributes(x)
    if (attr$missing) {
        list(type="H5T_NATIVE_INT8", placeholder=-1L, size=attr$non_zero)
    } else {
        list(type="H5T_NATIVE_INT8", placeholder=NULL, size=attr$non_zero)
    }
}
