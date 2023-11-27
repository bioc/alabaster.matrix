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
    range(unlist(lapply(collated, function(y) y[[name]])))
}

aggregate_any <- function(collated, name) {
    any(vapply(collated, function(y) y[[name]], TRUE))
}

aggregate_max <- function(collated, name) {
    max(unlist(lapply(collated, function(y) y[[name]])), na.rm=TRUE)
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

setMethod("collect_integer_attributes", "array", .simple_integer_collector)

setMethod("collect_integer_attributes", "ANY", function(x) {
    collated <- blockApply(x, .simple_integer_collector)
    list(
        range=aggregate_range(collated, "range"),
        missing=aggregate_any(collated, "missing")
    )
})

optimize_integer_storage <- function(x) {
    attr <- collect_integer_attributes(x)

    if (attr$missing) {
        lower <- attr$range[1]
        upper <- attr$range[2]
        if (is.infinite(lower)) {
            return(list(type="H5T_NATIVE_INT8", placeholder=as.integer(-2^7)))
        }

        if (lower < 0L) {
            if (lower > -2^7 && upper < 2^7) {
                return(list(type="H5T_NATIVE_INT8", placeholder=as.integer(-2^7)))
            } else if (lower > -2^15 && upper < 2^15) {
                return(list(type="H5T_NATIVE_INT16", placeholder=as.integer(-2^15)))
            }
        } else {
            if (upper < 2^8 - 1) {
                return(list(type="H5T_NATIVE_UINT8", placeholder=as.integer(2^8-1)))
            } else if (upper < 2^16 - 1) {
                return(list(type="H5T_NATIVE_UINT16", placeholder=as.integer(2^16-1)))
            }
        }

        return(list(type="H5T_NATIVE_INT32", placeholder=NA_integer_))

    } else {
        lower <- attr$range[1]
        upper <- attr$range[2]
        if (is.infinite(lower)) {
            return(list(type="H5T_NATIVE_INT8", placeholder=NULL))
        }

        if (lower < 0L) {
            if (lower >= -2^7 && upper < 2^7) {
                return(list(type="H5T_NATIVE_INT8", placeholder=NULL))
            } else if (lower >= -2^15 && upper < 2^15) {
                return(list(type="H5T_NATIVE_INT16", placeholder=NULL))
            }
        } else {
            if (upper < 2^8) {
                return(list(type="H5T_NATIVE_UINT8", placeholder=NULL))
            } else if (upper < 2^16) {
                return(list(type="H5T_NATIVE_UINT16", placeholder=NULL))
            }
        }

        return(list(type="H5T_NATIVE_INT32", placeholder=NULL))
    }
}

###################################################
###################################################

setGeneric("collect_float_attributes", function(x) standardGeneric("collect_float_attributes"))

setMethod("collect_float_attributes", "array", collect_double_attributes)

setMethod("collect_float_attributes", "ddenseMatrix", function(x) collect_double_attributes(x@x))

setMethod("collect_float_attributes", "ANY", function(x) {
    collated <- blockApply(x, collect_double_attributes)

    output <- list(range=aggregate_range(collated, "range"))
    for (n in c("missing", "non_integer", "has_NaN", "has_Inf", "has_nInf", "has_lowest", "has_highest")) {
        output[[n]] <- aggregate_any(collated, n)
    }

    output
})

optimize_float_storage <- function(x) {
    attr <- collect_float_attributes(x)

    if (attr$missing) {
        if (!attr$non_integer) {
            lower <- attr$range[1]
            upper <- attr$range[2]
            if (lower < 0L) {
                if (lower > -2^7 && upper < 2^7) {
                    return(list(type="H5T_NATIVE_INT8", placeholder=-2^7))
                } else if (lower > -2^15 && upper < 2^15) {
                    return(list(type="H5T_NATIVE_INT16", placeholder=-2^15))
                } else if (lower > -2^31 && upper < 2^31) {
                    return(list(type="H5T_NATIVE_INT32", placeholder=-2^31))
                }
            } else {
                if (upper < 2^8-1) {
                    return(list(type="H5T_NATIVE_UINT8", placeholder=2^8-1))
                } else if (upper < 2^16-1) {
                    return(list(type="H5T_NATIVE_UINT16", placeholder=2^16-1))
                } else if (upper < 2^32-1) {
                    return(list(type="H5T_NATIVE_UINT32", placeholder=2^32-1))
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
            u <- Reduce(union, blockApply(x, function(y) unique(as.vector(y))))
            placeholder <- chooseMissingPlaceholderForHdf5(u)
        }

        return(list(type="H5T_NATIVE_DOUBLE", placeholder=placeholder))

    } else {
        if (!attr$non_integer) {
            lower <- attr$range[1]
            upper <- attr$range[2]
            if (lower < 0L) {
                if (lower >= -2^7 && upper < 2^7) {
                    return(list(type="H5T_NATIVE_INT8", placeholder=NULL))
                } else if (lower >= -2^15 && upper < 2^15) {
                    return(list(type="H5T_NATIVE_INT16", placeholder=NULL))
                } else if (lower >= -2^31 && upper < 2^31) {
                    return(list(type="H5T_NATIVE_INT32", placeholder=NULL))
                }
            } else {
                if (upper < 2^8) {
                    return(list(type="H5T_NATIVE_UINT8", placeholder=NULL))
                } else if (upper < 2^16) {
                    return(list(type="H5T_NATIVE_UINT16", placeholder=NULL))
                } else if (upper < 2^32) {
                    return(list(type="H5T_NATIVE_UINT32", placeholder=NULL))
                }
            }
        }

        return(list(type="H5T_NATIVE_DOUBLE", placeholder=NULL))
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

optimize_boolean_storage <- function(x) {
    if (anyNA(x)) {
        list(type="H5T_NATIVE_INT8", placeholder=-1L)
    } else {
        list(type="H5T_NATIVE_INT8", placeholder=NULL)
    }
}
