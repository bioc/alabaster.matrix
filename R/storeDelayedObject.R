#' Store operations in a DelayedArray
#'
#' Store the delayed operations of a \linkS4class{DelayedArray} in a HDF5 file.
#'
#' @param x Any of the delayed operation classes from \pkg{DelayedArray}.
#' @param file String containing the path to a HDF5 file.
#' @param name String containing the name of the group to save into.
#' @param ... Arguments to be passed to specific methods.
#' @return The contents of \code{x} are saved to \code{file}, and \code{NULL} is invisibly returned.
#' 
#' @author Aaron Lun
#' @examples
#' library(DelayedArray)
#' X <- DelayedArray(matrix(runif(100), ncol=20))
#' Y <- cbind(X, DelayedArray::ConstantArray(value=50, c(5, 10)))
#'
#' library(rhdf5)
#' temp <- tempfile()
#' dir.create(temp)
#'
#' fpath <- file.path(temp, "foo.h5")
#' fhandle <- H5Fcreate(fpath)
#' storeDelayedObject(Y@seed, fhandle, "YAY")
#' rhdf5::h5ls(fhandle)
#' H5Fclose(fhandle)
#'
#' fhandle <- H5Fopen(fpath, "H5F_ACC_RDONLY")
#' reloadDelayedObject(fhandle, "YAY")
#' H5Fclose(fhandle)
#'
#' @aliases
#' storeDelayedObject
#' reloadDelayedObject
#' storeDelayedObject,ConstantArraySeed-method
#' storeDelayedObject,DelayedAbind-method
#' storeDelayedObject,ANY-method
#' storeDelayedObject,DelayedAperm-method
#' storeDelayedObject,DelayedNaryIsoOp-method
#' storeDelayedObject,DelayedSetDimnames-method
#' storeDelayedObject,DelayedSubassign-method
#' storeDelayedObject,DelayedSubset-method
#' storeDelayedObject,DelayedUnaryIsoOpStack-method
#' storeDelayedObject,DelayedUnaryIsoOpWithArgs-method
#' storeDelayedObject,SVT_SparseMatrix-method
#' storeDelayedObject,array-method
#' storeDelayedObject,denseMatrix-method
#' storeDelayedObject,sparseMatrix-method
#'
#' @name storeDelayedObject
NULL

r2value_type_mapping <- c(logical="BOOLEAN", integer="INTEGER", double="FLOAT", character="STRING")
to_value_type <- function(type) {
    if (!(type %in% names(r2value_type_mapping))) {
        stop("cannot map type '", type, "' to a value type")
    }
    r2value_type_mapping[[type]]
}

value2alabaster_type_mapping <- c(BOOLEAN="boolean", INTEGER="integer", FLOAT="number", STRING="string")
to_alabaster_type <- function(vtype) {
    if (!(vtype %in% names(value2alabaster_type_mapping))) {
        stop("cannot map type '", type, "' to an alabaster type")
    }
    value2alabaster_type_mapping[[vtype]]
}

#' @import alabaster.base rhdf5
load_vector_for_chihaya <- function(handle, name, version) {
    dhandle <- H5Dopen(handle, name)
    on.exit(H5Dclose(dhandle), add=TRUE, after=FALSE)
    contents <- H5Dread(dhandle, drop=TRUE) 
    if (is.raw(contents)) {
        storage.mode(contents) <- "integer"
    }

    type <- h5_read_attribute(dhandle, "type")
    type <- to_alabaster_type(type)
    missing.placeholder <- h5_read_attribute(dhandle, "missing_placeholder", check=TRUE, default=NULL)
    h5_cast(contents, expected.type=type, missing.placeholder=missing.placeholder)
}

#' @import alabaster.base rhdf5
save_vector_for_chihaya <- function(handle, name, x, version, scalar) {
    info <- transformVectorForHdf5(x)

    dtype <- NULL
    rtype <- typeof(x)
    if (rtype == "double") {
        dtype <- "H5T_NATIVE_DOUBLE"
    } else if (rtype == "integer") {
        dtype <- "H5T_NATIVE_INT32"
    } else if (rtype == "logical") {
        dtype <- "H5T_NATIVE_INT8"
    }

    dhandle <- h5_write_vector(handle, name, info$transformed, type=dtype, scalar=scalar, emit=TRUE)
    on.exit(H5Dclose(dhandle), add=TRUE, after=FALSE)

    if (!is.null(info$placeholder)) {
        h5_write_attribute(dhandle, "missing_placeholder", info$placeholder, type=dtype, scalar=TRUE)
    }
    h5_write_attribute(dhandle, "type", to_value_type(rtype), scalar=TRUE)

    invisible(NULL)
}

#######################################################
#######################################################

chihaya_array_registry <- list()
chihaya_operation_registry <- list()
chihaya_type_hint_registry <- list()

#' @export
reloadDelayedObject <- function(handle, name, version=package_version("1.1"), ...) {
    ghandle <- H5Gopen(handle, name)
    on.exit(H5Gclose(ghandle), add=TRUE, after=FALSE)

    output <- NULL
    if (h5_object_exists(ghandle, "_r_type_hint")) {
        hint <- h5_read_vector(ghandle, "_r_type_hint")
        if (hint %in% names(chihaya_type_hint_registry)) {
            FUN <- chihaya_type_hint_registry[[hint]]
            output <- tryCatch(
                FUN(ghandle, version=version, ...),
                error=function(e) NULL
            )
        }
    }

    if (is.null(output)) {
        objtype <- h5_read_attribute(ghandle, "delayed_type")
        if (objtype == "array") {
            arrtype <- h5_read_attribute(ghandle, "delayed_array")
            FUN <- chihaya_array_registry[[arrtype]]
        } else {
            optype <- h5_read_attribute(ghandle, "delayed_operation")
            FUN <- chihaya_operation_registry[[optype]]
        }

        if (is.null(FUN)) {
            stop(objtype, " type '", optype, "' is not yet supported")
        }
        output <- FUN(ghandle, version=version, ...)
    }

    if (!is(output, "DelayedArray")) {
        output <- DelayedArray(output)
    }
    output
}

#######################################################
#######################################################

#' @export
#' @import rhdf5
setMethod("storeDelayedObject", "ConstantArraySeed", function(x, handle, name, version=package_version("1.1"), ...) {
    ghandle <- H5Gcreate(handle, name)
    on.exit(H5Gclose(ghandle), add=TRUE, after=FALSE)

    h5_write_attribute(ghandle, "delayed_type", "array", scalar=TRUE)
    h5_write_attribute(ghandle, "delayed_array", "constant array", scalar=TRUE)
    h5_write_vector(ghandle, "dimensions", dim(x), compress=0, type="H5T_NATIVE_UINT32")
    save_vector_for_chihaya(ghandle, "value", x@value, version=version, scalar=TRUE)
})

#' @import rhdf5 DelayedArray
chihaya_array_registry[["constant array"]] <- function(handle, version, ...) {
    dim <- h5_read_vector(handle, "dimensions")
    val <- load_vector_for_chihaya(handle, "value", version=version)
    ConstantArray(dim, value=val)
}

#######################################################
#######################################################

#' @import rhdf5 alabaster.base
save_dense_array_for_chihaya <- function(x, handle, name, extract.native, version=package_version("1.1"), ...) {
    ghandle <- H5Gcreate(handle, name)
    on.exit(H5Gclose(ghandle), add=TRUE, after=FALSE)

    h5_write_attribute(ghandle, "delayed_type", "array", scalar=TRUE)
    h5_write_attribute(ghandle, "delayed_array", "dense array", scalar=TRUE)

    optimized <- optimize_storage(x)
    h5_write_array(
        ghandle, 
        name="data", 
        x=x, 
        type=optimized$type, 
        placeholder=optimized$placeholder, 
        extract.native=extract.native
    )

    dhandle <- H5Dopen(ghandle, "data")
    on.exit(H5Dclose(dhandle), add=TRUE, after=FALSE)
    h5_write_attribute(dhandle, "type", to_value_type(type(x)), scalar=TRUE)
    if (!is.null(optimized$placeholder)) {
        h5_write_attribute(dhandle, "missing_placeholder", optimized$placeholder, type=optimized$type, scalar=TRUE)
    }

    h5_write_vector(ghandle, "native", 0L, type="H5T_NATIVE_INT8", scalar=TRUE)

    dn <- dimnames(x)
    if (!is.null(dn)) {
        save_dimnames_for_chihaya(ghandle, rev(dn))
    }

    invisible(NULL)
}

#' @export
setMethod("storeDelayedObject", "array", function(x, handle, name, version=package_version("1.1"), save.external.array=FALSE, ...) {
    if (save.external.array) {
        return(callNextMethod()) # calls the ANY method
    }
    save_dense_array_for_chihaya(x, handle, name, extract.native=identity, version=version, ...)
})

#' @export
setMethod("storeDelayedObject", "denseMatrix", function(x, handle, name, version=package_version("1.1"), save.external.array=FALSE, ...) {
    if (save.external.array) {
        return(callNextMethod()) # calls the ANY method.
    } 

    extract.native <- NULL
    if (is(x, "dMatrix") || is(x, "lMatrix")) {
        extract.native <- function(y) y@x
    }
    save_dense_array_for_chihaya(x, handle, name, extract.native=extract.native, version=version, ...)
})

#' @import rhdf5 DelayedArray
chihaya_array_registry[["dense array"]] <- function(handle, version, ...) {
    dhandle <- H5Dopen(handle, "data")
    on.exit(H5Dclose(dhandle), add=TRUE, after=FALSE)
    data <- H5Dread(dhandle)
    if (is.raw(data)) {
        storage.mode(data) <- "integer"
    }

    type <- h5_read_attribute(dhandle, "type")
    type <- to_alabaster_type(type)
    missing.placeholder <- h5_read_attribute(dhandle, "missing_placeholder", check=TRUE, default=NULL)
    data <- h5_cast(data, expected.type=type, missing.placeholder=missing.placeholder)

    if (h5_object_exists(handle, "dimnames")) {
        dn <- load_dimnames_for_chihaya(handle)
        dimnames(data) <- rev(dn)
    }

    if (h5_read_vector(handle, "native") == 1L) {
        data <- t(data)
    }

    data
}

#######################################################
#######################################################

#' @import rhdf5
save_sparse_matrix_for_chihaya <- function(x, handle, name, version=package_version("1.1"), ...) {
    ghandle <- H5Gcreate(handle, name)
    on.exit(H5Gclose(ghandle), add=TRUE, after=FALSE)

    h5_write_attribute(ghandle, "delayed_type", "array", scalar=TRUE)
    h5_write_attribute(ghandle, "delayed_array", "sparse matrix", scalar=TRUE)

    h5_write_vector(ghandle, "shape", dim(x), type="H5T_NATIVE_UINT32")

    optimized <- optimize_storage(x)
    h5_write_sparse_matrix(x, handle=ghandle, details=optimized)

    dhandle <- H5Dopen(ghandle, "data")
    on.exit(H5Dclose(dhandle), add=TRUE, after=FALSE)
    h5_write_attribute(dhandle, "type", to_value_type(type(x)), scalar=TRUE)
    if (!is.null(optimized$placeholder)) {
        h5_write_attribute(dhandle, missingPlaceholderName, optimized$placeholder, type=optimized$type, scalar=TRUE)
    }

    # Better if we didn't save the layout at all, but whatever.
    col <- h5_read_attribute(ghandle, "layout") == "CSC"
    h5_write_vector(ghandle, "by_column", as.integer(col), type="H5T_NATIVE_INT8", scalar=TRUE)
    H5Adelete(ghandle, "layout") 

    dn <- dimnames(x)
    if (!is.null(dn)) {
        save_dimnames_for_chihaya(ghandle, dn)
    }

    invisible(NULL)
}

#' @export
setMethod("storeDelayedObject", "sparseMatrix", function(x, handle, name, version=package_version("1.1"), save.external.array=FALSE, ...) {
    if (save.external.array) {
        return(callNextMethod()) # calls the ANY method
    }
    save_sparse_matrix_for_chihaya(x, handle, name, version=version, ...)
})

#' @export
setMethod("storeDelayedObject", "SVT_SparseMatrix", function(x, handle, name, version=package_version("1.1"), save.external.array=FALSE, ...) {
    if (save.external.array) {
        return(callNextMethod()) # calls the ANY method
    }
    save_sparse_matrix_for_chihaya(x, handle, name, version=version, ...)
})

#' @import rhdf5 DelayedArray
chihaya_array_registry[["sparse matrix"]] <- function(handle, version, ...) {
    indices <- h5_read_vector(handle, "indices")
    indptr <- h5_read_vector(handle, "indptr")
    data <- load_vector_for_chihaya(handle, "data", version=version)

    svt <- vector("list", length(indptr) - 1L)
    for (i in seq_along(svt)) {
        idx <- indptr[i] + seq_len(indptr[i+1] - indptr[i])
        if (length(idx)) {
            # As of version >= 1, data is first and indices are second.
            svt[[i]] <- list(data[idx], indices[idx])
        }
    }

    dim <- h5_read_vector(handle, "shape")
    dn <- list(NULL, NULL)
    if (h5_object_exists(handle, "dimnames")) {
        dn <- load_dimnames_for_chihaya(handle)
    }

    transposed <- (h5_read_vector(handle, "by_column") == 0L)
    if (transposed) {
        dn <- rev(dn)
        dim <- rev(dim)
    }

    output <- new("SVT_SparseMatrix", SVT=svt, dim=dim, type=typeof(data), dimnames=dn)
    if (transposed) {
        output <- t(output)
    }

    output
}

#######################################################
#######################################################

#' @export
#' @import rhdf5
setMethod("storeDelayedObject", "DelayedAbind", function(x, handle, name, version=package_version("1.1"), ...) {
    ghandle <- H5Gcreate(handle, name)
    on.exit(H5Gclose(ghandle), add=TRUE, after=FALSE)

    h5_write_attribute(ghandle, "delayed_type", "operation", scalar=TRUE)
    h5_write_attribute(ghandle, "delayed_operation", "combine", scalar=TRUE)
    h5_write_vector(ghandle, "along", x@along - 1L, type="H5T_NATIVE_UINT32", scalar=TRUE)

    shandle <- H5Gcreate(ghandle, "seeds")
    on.exit(H5Gclose(shandle), add=TRUE, after=FALSE)
    h5_write_attribute(shandle, "length", length(x@seeds), type="H5T_NATIVE_UINT32", scalar=TRUE)

    for (i in seq_along(x@seeds)) {
        storeDelayedObject(x@seeds[[i]], shandle, as.character(i - 1L), version=version, ...)
    }

    invisible(NULL)
})

#' @import alabaster.base rhdf5 DelayedArray
chihaya_operation_registry[["combine"]] <- function(handle, version, ...) {
    shandle <- H5Gopen(handle, "seeds")
    on.exit(H5Gclose(shandle), add=TRUE, after=FALSE)

    len <- h5_read_attribute(shandle, "length")
    along <- h5_read_vector(handle, "along")

    seeds <- vector("list", len)
    for (i in seq_len(len)) {
        seeds[[i]] <- reloadDelayedObject(shandle, as.character(i - 1L), version=version, ...)
    }

    if (along == 0L) {
        do.call(arbind, seeds)
    } else {
        do.call(acbind, seeds)
    }
}

#######################################################
#######################################################

#' @export
#' @import rhdf5
setMethod("storeDelayedObject", "DelayedAperm", function(x, handle, name, version=package_version("1.1"), ...) {
    ghandle <- H5Gcreate(handle, name)
    on.exit(H5Gclose(ghandle), add=TRUE, after=FALSE)

    h5_write_attribute(ghandle, "delayed_type", "operation", scalar=TRUE)
    h5_write_attribute(ghandle, "delayed_operation", "transpose", scalar=TRUE)

    h5_write_vector(ghandle, "permutation", x@perm - 1L, type="H5T_NATIVE_UINT32")
    storeDelayedObject(x@seed, ghandle, "seed", version=version, ...)
    invisible(NULL)
})

#' @import DelayedArray rhdf5
chihaya_operation_registry[["transpose"]] <- function(handle, version, ...) {
    x <- reloadDelayedObject(handle, "seed", version=version, ...)
    perm <- h5_read_vector(handle, "permutation")
    aperm(x, perm + 1L)
}

#######################################################
#######################################################

chihaya.supported.Arith <- c("+", "-", "*", "^", "/", "%%", "%/%")
chihaya.supported.Compare <- c("==", ">", "<", "!=", "<=", ">=")
chihaya.supported.Logic <- c("&", "|")
chihaya.supported.Ops <- c(chihaya.supported.Arith, chihaya.supported.Compare, chihaya.supported.Logic)

translate_logic_Ops_for_chihaya <- function(method) {
    if (method == "&") {
        "&&"
    } else if (method == "|") {
        "||"
    } else {
        method
    }
}

translate_logic_Ops_from_chihaya <- function(method) {
    if (method == "&&") {
        "&"
    } else if (method == "||") {
        "|"
    } else {
        method
    }
}

#' @export
#' @import rhdf5
setMethod("storeDelayedObject", "DelayedNaryIsoOp", function(x, handle, name, version=package_version("1.1"), ...) {
    ghandle <- H5Gcreate(handle, name)
    on.exit(H5Gclose(ghandle), add=TRUE, after=FALSE)
    h5_write_attribute(ghandle, "delayed_type", "operation", scalar=TRUE)

    # Figuring out the identity of the operation.
    chosen <- NULL
    for (p in chihaya.supported.Ops) {
        if (identical(x@OP, get(p, envir=baseenv()))) {
            chosen <- p
            break
        }
    }
    if (is.null(chosen)) {
        stop("unknown operation in ", class(x))
    }

    if (chosen %in% chihaya.supported.Arith) {
        op <- "binary arithmetic"
    } else if (chosen %in% chihaya.supported.Compare) {
        op <- "binary comparison"
    } else if (chosen %in% chihaya.supported.Logic) {
        op <- "binary logic"
        chosen <- translate_logic_Ops_for_chihaya(chosen)
    }
    h5_write_attribute(ghandle, "delayed_operation", op, scalar=TRUE)
    h5_write_vector(ghandle, "method", chosen, scalar=TRUE)

    if (length(x@seeds) != 2) {
        stop("expected exactly two seeds for 'DelayedNaryIsoOp'")
    }
    if (length(x@Rargs)) {
        stop("expected no additional right arguments for 'DelayedNaryIsoOp'")
    }

    storeDelayedObject(x@seeds[[1]], ghandle, "left")
    storeDelayedObject(x@seeds[[2]], ghandle, "right")
    invisible(NULL)
})

#' @import DelayedArray
chihaya_load_binary_op <- function(handle, version, logic, ...) {
    left <- reloadDelayedObject(handle, "left", version=version, ...)
    right <- reloadDelayedObject(handle, "right", version=version, ...)
    op <- h5_read_vector(handle, "method")
    if (logic) {
        op <- translate_logic_Ops_from_chihaya(op)
    }
    get(op, envir=baseenv())(left, right)
}

chihaya_operation_registry[["binary arithmetic"]] <- function(handle, version, ...) chihaya_load_binary_op(handle, version, logic=FALSE, ...)
chihaya_operation_registry[["binary comparison"]] <- function(handle, version, ...) chihaya_load_binary_op(handle, version, logic=FALSE, ...)
chihaya_operation_registry[["binary logic"]] <- function(handle, version, ...) chihaya_load_binary_op(handle=handle, version=version, logic=TRUE, ...)

#######################################################
#######################################################

#' @import rhdf5 alabaster.base
save_dimnames_for_chihaya <- function(handle, dimnames) {
    dhandle <- H5Gcreate(handle, "dimnames")
    on.exit(H5Gclose(dhandle), add=TRUE, after=FALSE)
    h5_write_attribute(dhandle, "length", length(dimnames), type="H5T_NATIVE_UINT32", scalar=TRUE)

    for (i in seq_along(dimnames)) {
        dn <- dimnames[[i]]
        if (is.character(dn)) { # avoid NULLs, -1's.
            h5_write_vector(dhandle, as.character(i - 1L), dn)
        }
    }
}

#' @import rhdf5 alabaster.base
load_dimnames_for_chihaya <- function(handle) {
    dhandle <- H5Gopen(handle, "dimnames")
    on.exit(H5Gclose(dhandle), add=TRUE, after=FALSE)

    dlen <- h5_read_attribute(dhandle, "length")
    dnames <- vector("list", dlen)
    for (i in seq_along(dnames)) {
        n <- as.character(i - 1L)
        if (h5_object_exists(dhandle, n)) {
            dnames[[i]] <- h5_read_vector(dhandle, n)
        }
    }

    dnames
}

#' @export
#' @import rhdf5
setMethod("storeDelayedObject", "DelayedSetDimnames", function(x, handle, name, version=package_version('1.1'), ...) {
    ghandle <- H5Gcreate(handle, name)
    on.exit(H5Gclose(ghandle), add=TRUE, after=FALSE)

    h5_write_attribute(ghandle, "delayed_type", "operation", scalar=TRUE)
    h5_write_attribute(ghandle, "delayed_operation", "dimnames", scalar=TRUE)
    save_dimnames_for_chihaya(ghandle, x@dimnames)

    storeDelayedObject(x@seed, ghandle, "seed", version=version, ...)
    invisible(NULL)
})

#' @import rhdf5 DelayedArray
chihaya_operation_registry[["dimnames"]] <- function(handle, version, ...) {
    x <- reloadDelayedObject(handle, "seed", version=version, ...)
    dimnames(x) <- load_dimnames_for_chihaya(handle)
    x
}

#######################################################
#######################################################

#' @import rhdf5 alabaster.base
save_chihaya_indices <- function(handle, name, indices) { 
    ihandle <- H5Gcreate(handle, name)
    on.exit(H5Gclose(ihandle), add=TRUE, after=FALSE)
    h5_write_attribute(ihandle, "length", length(indices), type="H5T_NATIVE_UINT32", scalar=TRUE)

    for (i in seq_along(indices)) {
        ii <- indices[[i]]
        if (!is.null(ii)) {
            h5_write_vector(ihandle, as.character(i - 1L), ii - 1L, type="H5T_NATIVE_UINT32")
        }
    }
}

#' @import rhdf5 alabaster.base
load_chihaya_indices <- function(handle, name) {
    ihandle <- H5Gopen(handle, name)
    on.exit(H5Gclose(ihandle), add=TRUE, after=FALSE)
    ilen <- h5_read_attribute(ihandle, "length")

    indices <- vector("list", ilen)
    for (i in seq_len(ilen)) {
        n <- as.character(i - 1L)
        if (h5_object_exists(ihandle, n)) {
            indices[[i]] <- h5_read_vector(ihandle, n) + 1L
        } else {
            indices[[i]] <- substitute()
        }
    }

    indices
}

#' @export
#' @import rhdf5
setMethod("storeDelayedObject", "DelayedSubassign", function(x, handle, name, version=package_version("1.1"), ...) {
    ghandle <- H5Gcreate(handle, name)
    on.exit(H5Gclose(ghandle), add=TRUE, after=FALSE)

    h5_write_attribute(ghandle, "delayed_type", "operation", scalar=TRUE)
    h5_write_attribute(ghandle, "delayed_operation", "subset assignment", scalar=TRUE)

    save_chihaya_indices(ghandle, "index", x@Lindex) 
    storeDelayedObject(x@seed, ghandle, "seed", version=version, ...)
    storeDelayedObject(x@Rvalue, ghandle, "value", version=version, ...)
    invisible(NULL)
})

#' @import rhdf5 DelayedArray
chihaya_operation_registry[["subset assignment"]] <- function(handle, version, ...) {
    x <- reloadDelayedObject(handle, "seed", version=version, ...)
    value <- reloadDelayedObject(handle, "value", version=version, ...)
    indices <- load_chihaya_indices(handle, "index")
    do.call(`[<-`, c(list(x=x), indices, list(value=value)))
} 

#######################################################
#######################################################

#' @export
#' @import rhdf5
setMethod("storeDelayedObject", "DelayedSubset", function(x, handle, name, version=package_version("1.1"), ...) {
    ghandle <- H5Gcreate(handle, name)
    on.exit(H5Gclose(ghandle), add=TRUE, after=FALSE)

    h5_write_attribute(ghandle, "delayed_type", "operation", scalar=TRUE)
    h5_write_attribute(ghandle, "delayed_operation", "subset", scalar=TRUE)

    save_chihaya_indices(ghandle, "index", x@index) 
    storeDelayedObject(x@seed, ghandle, "seed", version=version, ...)
    invisible(NULL)
})

#' @import rhdf5 DelayedArray
chihaya_operation_registry[["subset"]] <- function(handle, version, ...) {
    x <- reloadDelayedObject(handle, "seed", version=version, ...)
    indices <- load_chihaya_indices(handle, "index")
    do.call(`[`, c(list(x), indices, list(drop=FALSE)))
} 

#######################################################
#######################################################

chihaya.dump.unary.Math <- function(handle, OP) {
    envir <- environment(OP)
    generic <- envir$`.Generic`

    if (!is.null(generic)) {
        direct.support <- c("abs", "sign", "sqrt", "ceiling", "floor", "trunc",
                            "exp", "expm1", "log1p", "cos", "cosh", "sin",
                            "sinh", "tan", "tanh", "acos", "acosh", "asin",
                            "asinh", "atan", "atanh", "gamma", "lgamma",
                            "digamma", "trigamma")

        if (generic %in% direct.support) {
            h5_write_attribute(handle, "delayed_operation", 'unary math', scalar=TRUE)
            h5_write_vector(handle, "method", generic, scalar=TRUE)
            return(TRUE)
        }

        log.base.support <- c(log2=2, log10=10)
        if (generic %in% names(log.base.support)) {
            h5_write_attribute(handle, "delayed_operation", 'unary math', scalar=TRUE)
            h5_write_vector(handle, "method", "log", scalar=TRUE)
            h5_write_vector(handle, "base", log.base.support[[generic]], type="H5T_NATIVE_DOUBLE", scalar=TRUE)
            return(TRUE)
        }
    }

    # Special case for the general case log.
    if (isTRUE(all.equal(as.character(body(OP)), c("log", "a", "base")))) {
        h5_write_attribute(handle, "delayed_operation", 'unary math', scalar=TRUE)
        h5_write_vector(handle, "method", "log", scalar=TRUE)

        base <- envir$base
        if (base != exp(1)) {
            h5_write_vector(handle, "base", base, type="H5T_NATIVE_DOUBLE", scalar=TRUE)
        }
        return(TRUE)
    }

    FALSE
}

chihaya.dump.unary.Math2 <- function(handle, OP) {
    envir <- environment(OP)
    generic <- envir$`.Generic`

    if (generic %in% getGroupMembers("Math2")) {
        h5_write_attribute(handle, "delayed_operation", 'unary math', scalar=TRUE)
        h5_write_vector(handle, "method", generic, scalar=TRUE)
        h5_write_vector(handle, "digits", envir$digits, type="H5T_NATIVE_INT32", scalar=TRUE)
        return(TRUE)
    }

    FALSE
}

chihaya.dump.unary.Ops <- function(handle, OP) {
    envir <- environment(OP)
    generic <- envir$`.Generic`

    if (!(generic %in% chihaya.supported.Ops)) {
        return(FALSE)
    }

    delayed_op <- NULL
    if (generic %in% chihaya.supported.Arith) {
        delayed_op <- "unary arithmetic"
    } else if (generic %in% chihaya.supported.Compare) {
        delayed_op <- "unary comparison"
    } else if (generic %in% chihaya.supported.Logic) {
        delayed_op <- "unary logic"
        generic <- translate_logic_Ops_for_chihaya(generic)
    }
    h5_write_attribute(handle, "delayed_operation", delayed_op, scalar=TRUE)
    h5_write_vector(handle, "method", generic, scalar=TRUE)

    e1 <- envir$e1
    e2 <- envir$e2

    if (missing(e2)) {
        if (!generic %in% c("+", "-")) {
            stop("second argument can only be missing for unary '+' or '-'")
        }
        h5_write_vector(handle, "side", "none", scalar=TRUE)

    } else {
        right <- is(e1, "DelayedArray") # i.e., is the operation applied to the left of the seed?
        left <- is(e2, "DelayedArray") # i.e., is the operation applied to the left of the seed?

        h5_write_vector(handle, "side", if (left) "left" else "right", scalar=TRUE)
        val <- if (left) e1 else e2
        if (length(val) == 1) {
            save_vector_for_chihaya(handle, "value", val, version=version, scalar=TRUE)
        } else {
            # Don't think this ever gets called, because otherwise
            # we'd be dealing with a DelayedIsoOpWithArgs.
            # Nonetheless, we'll throw in the necessary code.
            h5_write_vector(handle, "along", 0, type="H5T_NATIVE_UINT32", scalar=TRUE)
            save_vector_for_chihaya(handle, "value", val, version=version)
        }
    }

    TRUE
}

chihaya.unary.logic.Ops <- c(is.infinite="is_infinite", is.nan="is_nan", is.finite="is_finite")

chihaya.dump.unary.other <- function(handle, OP) {
    envir <- environment(OP)
    generic <- envir$`.Generic`

    if (generic == "!") {
        h5_write_attribute(handle, "delayed_operation", "unary logic", scalar=TRUE)
        h5_write_vector(handle, "method", "!", scalar=TRUE)
        return(TRUE)

    } else if (generic %in% names(chihaya.unary.logic.Ops)) {
        h5_write_attribute(handle, "delayed_operation", "unary special check", scalar=TRUE)
        h5_write_vector(handle, "method", chihaya.unary.logic.Ops[[generic]], scalar=TRUE)
        return(TRUE)
    }

    FALSE
}

save_iso_op_stack <- function(handle, name, OPS, i, seed, version, ...) {
    if (i == 0L) {
        storeDelayedObject(seed, handle, name, version=version, ...)
        return(NULL)
    }

    OP <- OPS[[i]]
    ghandle <- H5Gcreate(handle, name)
    on.exit(H5Gclose(ghandle), add=TRUE, after=FALSE)
    h5_write_attribute(ghandle, "delayed_type", 'operation', scalar=TRUE)

    if (chihaya.dump.unary.Math(ghandle, OP)) {
        return(save_iso_op_stack(ghandle, "seed", OPS, i - 1L, seed=seed, version=version, ...))
    }

    if (chihaya.dump.unary.Math2(ghandle, OP)) {
        return(save_iso_op_stack(ghandle, "seed", OPS, i - 1L, seed=seed, version=version, ...))
    }

    if (chihaya.dump.unary.Ops(ghandle, OP)) {
        return(save_iso_op_stack(ghandle, "seed", OPS, i - 1L, seed=seed, version=version, ...))
    }

    if (chihaya.dump.unary.other(ghandle, OP)) {
        return(save_iso_op_stack(ghandle, "seed", OPS, i - 1L, seed=seed, version=version, ...))
    }

    stop("unknown OPS[[", i, "]] function when saving a 'DelayedUnaryIsoOpStack'")
}

#' @export
setMethod("storeDelayedObject", "DelayedUnaryIsoOpStack", function(x, handle, name, version=package_version("1.1"), ...) {
    # This saves in reverse order, as first operation is first applied (and thus needs to be closer to the leaf of the delayed tree).
    save_iso_op_stack(handle, name, OPS=x@OPS, i=length(x@OPS), seed=x@seed, version=version, ...)
    invisible(NULL)
})

#' @export
setMethod("storeDelayedObject", "DelayedUnaryIsoOpWithArgs", function(x, handle, name, version=package_version("1.1"), ...) {
    ghandle <- H5Gcreate(handle, name)
    on.exit(H5Gclose(ghandle), add=TRUE, after=FALSE)
    h5_write_attribute(ghandle, "delayed_type", 'operation', scalar=TRUE)

    # Figuring out the identity of the operation.
    chosen <- NULL
    for (p in chihaya.supported.Ops) {
        if (identical(x@OP, get(p, envir=baseenv()))) {
            chosen <- p
            break
        }
    }
    if (is.null(chosen)) {
        stop("unknown operation in ", class(x))
    }

    op <- NULL
    if (chosen %in% chihaya.supported.Logic) {
        op <- "unary logic"
        chosen <- translate_logic_Ops_for_chihaya(chosen)
    } else if (chosen %in% chihaya.supported.Compare) {
        op <- "unary comparison"
    } else if (chosen %in% chihaya.supported.Arith) {
        op <- "unary arithmetic"
    }
    h5_write_attribute(ghandle, "delayed_operation", op, scalar=TRUE)
    h5_write_vector(ghandle, "method", chosen, scalar=TRUE)

    # Saving the left and right args. There should only be one or the other.
    # as the presence of both is not commutative.
    if (length(x@Rargs) + length(x@Largs) !=1) {
        stop("'DelayedUnaryIsoApWithArgs' should operate on exactly one argument")
    }

    left <- length(x@Largs) > 0
    if (left) {
        args <- x@Largs[[1]]
        along <- x@Lalong[1]
    } else {
        args <- x@Rargs[[1]]
        along <- x@Ralong[1]
    }

    h5_write_vector(ghandle, "side", if (left) "left" else "right", scalar=TRUE)
    h5_write_vector(ghandle, "along", along - 1L, type="H5T_NATIVE_UINT32", scalar=TRUE)

    if (length(args) == 1L) {
        save_vector_for_chihaya(ghandle, "value", args, version=version, scalar=TRUE)
    } else if (is.null(dim(args)) || length(dim(args)) == 1L) {
        save_vector_for_chihaya(ghandle, "value", args, version=version)
    } else {
        stop("multi-dimensional 'value' not supported in 'DelayedUnaryIsoOpWithArgs'")
    }

    storeDelayedObject(x@seed, ghandle, "seed", version=version, ...) 
    invisible(NULL)
})

chihaya_operation_registry[["unary math"]] <- function(handle, version, ...) {
    x <- reloadDelayedObject(handle, "seed", version=version, ...)
    method <- h5_read_vector(handle, "method")
    output <- NULL

    if (method == "log") {
        if (h5_object_exists(handle, "base")){
            base <- h5_read_vector(handle, "base")
            output <- log(x, base=base)
        } else {
            output <- log(x)
        }

    } else if (method %in% getGroupMembers("Math2")) {
        digits <- h5_read_vector(handle, "digits")
        FUN <- get(method, envir=baseenv())
        output <- FUN(x, digits=digits)

    } else {
        FUN <- get(method, envir=baseenv())
        output <- FUN(x)
    }

    output
}

apply_unary_op_with_value <- function(x, op, side, handle, version) {
    value <- load_vector_for_chihaya(handle, "value", version=version)
    FUN <- get(op, envir=baseenv())

    simple <- FALSE
    if (length(value) == 1L) {
        simple <- TRUE
    } else {
        along <- h5_read_vector(handle, "along")
        simple <- along == 0
    }

    output <- NULL
    if (simple) {
        if (side == "left") {
            output <- FUN(value, x)
        } else if (side == "right") {
            output <- FUN(x, value)
        } else {
            stop("unknown side '", side, "'")
        }

    } else {
        # Stolen from base::sweep.
        along <- along + 1L
        perm <- c(along, seq_along(dim(x))[-along])
        tmp <- aperm(x, perm)

        if (side == "left") {
            tmp <- FUN(value, tmp)
        } else if (side == "right") {
            tmp <- FUN(tmp, value)
        } else {
            stop("unknown side '", side, "'")
        }

        output <- aperm(tmp, order(perm))
    }

    output
}

chihaya_operation_registry[["unary logic"]] <- function(handle, version, ...) {
    x <- reloadDelayedObject(handle, "seed", version=version, ...)
    method <- h5_read_vector(handle, "method")

    output <- NULL
    if (method == "!") {
        output <- !x
    } else {
        method <- translate_logic_Ops_from_chihaya(method)
        side <- h5_read_vector(handle, "side")
        output <- apply_unary_op_with_value(x, op=method, side=side, handle=handle, version=version)
    }

    output
}

chihaya_operation_registry[["unary comparison"]] <- function(handle, version, ...) {
    x <- reloadDelayedObject(handle, "seed", version=version, ...)
    method <- h5_read_vector(handle, "method")
    side <- h5_read_vector(handle, "side")
    apply_unary_op_with_value(x, op=method, side=side, handle=handle, version=version)
}

chihaya_operation_registry[["unary special check"]] <- function(handle, version, ...) {
    x <- reloadDelayedObject(handle, "seed", version=version, ...)
    method <- h5_read_vector(handle, "method")
    method <- sub("_", ".", method)
    get(method, envir=baseenv())(x)
}

chihaya_operation_registry[["unary arithmetic"]] <- function(handle, version, ...) {
    x <- reloadDelayedObject(handle, "seed", version=version, ...)
    method <- h5_read_vector(handle, "method")
    side <- h5_read_vector(handle, "side")

    output <- NULL
    if (side == "none") {
        if (method == "+") {
            output <- x
        } else if (method == "-") {
            output <- -x
        } else {
            stop("unsupported method '", method, "' with side 'none'")
        }
    } else {
        output <- apply_unary_op_with_value(x, op=method, side=side, handle=handle, version=version)
    }

    output
}

#######################################################
#######################################################

#' @export
#' @import rhdf5
setMethod("storeDelayedObject", "ANY", function(x, handle, name, version=package_version("1.1"), ...) {
    ghandle <- H5Gcreate(handle, name)
    on.exit(H5Gclose(ghandle), add=TRUE, after=FALSE)

    if (is(x, "LowRankMatrixSeed")) { # From BiocSingular.
        h5_write_attribute(ghandle, "delayed_type", "operation", scalar=TRUE)
        h5_write_attribute(ghandle, "delayed_operation", "matrix product", scalar=TRUE)

        storeDelayedObject(x@rotation, ghandle, "left_seed", version=version, ...)
        h5_write_vector(ghandle, "left_orientation", "N", scalar=TRUE)
        storeDelayedObject(x@components, ghandle, "right_seed", version=version, ...)
        h5_write_vector(ghandle, "right_orientation", "T", scalar=TRUE)

    } else if (is(x, "ResidualMatrixSeed")) {
        h5_write_attribute(ghandle, "delayed_type", "operation", scalar=TRUE)
        h5_write_vector(ghandle, "_r_type_hint", "residual matrix", scalar=TRUE)

        # Mimic a transposition operation.
        xhandle <- ghandle
        if (x@transposed) {
            h5_write_attribute(ghandle, "delayed_operation", "transpose", scalar=TRUE)
            h5_write_vector(ghandle, "permutation", c(1L, 0L), type="H5T_NATIVE_UINT32") 
            xhandle <- H5Gcreate(ghandle, "seed")
            on.exit(H5Gclose(xhandle), add=TRUE, after=FALSE)
            h5_write_attribute(xhandle, "delayed_type", "operation", scalar=TRUE)
        }

        # Mimic a binary subtraction.
        h5_write_attribute(xhandle, "delayed_operation", "binary arithmetic", scalar=TRUE)
        h5_write_vector(xhandle, "method", "-", scalar=TRUE)
        storeDelayedObject(x@.matrix, xhandle, "left", version=version, ...)

        # Mimic a matrix product.
        rhandle <- H5Gcreate(xhandle, "right")
        on.exit(H5Gclose(rhandle), add=TRUE, after=FALSE)
        h5_write_attribute(rhandle, "delayed_type", "operation", scalar=TRUE)
        h5_write_attribute(rhandle, "delayed_operation", "matrix product", scalar=TRUE)
        storeDelayedObject(x@Q, rhandle, "left_seed", version=version, ...)
        h5_write_vector(rhandle, "left_orientation", "N", scalar=TRUE)
        storeDelayedObject(x@Qty, rhandle, "right_seed", version=version, ...)
        h5_write_vector(rhandle, "right_orientation", "N", scalar=TRUE)

    } else {
        h5_write_attribute(ghandle, "delayed_type", "array", scalar=TRUE)
        h5_write_attribute(ghandle, "delayed_array", "custom takane seed array", scalar=TRUE)

        exdir <- file.path(dirname(H5Fget_name(handle)), "seeds")
        dir.create(exdir, showWarnings=FALSE)
        n <- length(list.files(exdir))
        saveObject(x, file.path(exdir, n), ...)

        h5_write_vector(ghandle, "dimensions", dim(x), type="H5T_NATIVE_UINT32", compress=0)
        h5_write_vector(ghandle, "type", to_value_type(type(x)), scalar=TRUE)
        h5_write_vector(ghandle, "index", n, type="H5T_NATIVE_UINT32", scalar=TRUE)
    }
})

#' @import alabaster.base rhdf5 DelayedArray
chihaya_array_registry[["custom takane seed array"]] <- function(handle, version, custom.takane.realize=FALSE, ...) {
    index <- h5_read_vector(handle, "index")
    out <- readObject(file.path(dirname(H5Fget_name(handle)), "seeds", index), ...)

    if (custom.takane.realize) {
        if (is_sparse(out)) {
            out <- as(out, "SVT_SparseArray")
        } else {
            out <- as.array(out)
        }
    }

    out
}

#' @importFrom Matrix t
chihaya_operation_registry[["matrix product"]] <- function(handle, version, ...) {
    L <- as.matrix(reloadDelayedObject(handle, "left_seed"))
    Lori <- h5_read_vector(handle, "left_orientation")
    if (length(Lori) == 1 && as.character(Lori) == "T") {
        L <- t(L)
    } 

    R <- as.matrix(reloadDelayedObject(handle, "right_seed"))
    Rori <- h5_read_vector(handle, "right_orientation")
    if (length(Rori) == 1 && as.character(Rori) == "N") {
        R <- t(R)
    } 

    BiocSingular::LowRankMatrix(L, R)
}

chihaya_type_hint_registry[["residual matrix"]] <- function(handle, version, ...) {
    if (!isNamespaceLoaded("ResidualMatrix")) {
        loadNamespace("ResidualMatrix")
    }

    stopifnot(identical(h5_read_attribute(handle, "delayed_type"), "operation"))
    optype <- h5_read_attribute(handle, "delayed_operation")

    transposed <- FALSE
    xhandle <- handle
    if (optype == "transpose") {
        transposed <- TRUE
        xhandle <- H5Gopen(handle, "seed")
        on.exit(H5Gclose(xhandle), add=TRUE, after=FALSE)
        stopifnot(identical(h5_read_attribute(xhandle, "delayed_type"), "operation"))
        optype <- h5_read_attribute(xhandle, "delayed_operation")
    }

    stopifnot(identical(optype, "binary arithmetic"))
    stopifnot(identical(h5_read_vector(xhandle, "method"), "-"))
    .matrix <- reloadDelayedObject(xhandle, "left")
    .matrix <- .matrix@seed # can't be a DelayedArray inside ResidualMatrix, for various reasons...

    rhandle <- H5Gopen(xhandle, "right")
    on.exit(H5Gclose(rhandle), add=TRUE, after=FALSE)
    stopifnot(identical(h5_read_attribute(rhandle, "delayed_type"), "operation"))
    stopifnot(identical(h5_read_attribute(rhandle, "delayed_operation"), "matrix product"))
    stopifnot(identical(h5_read_vector(rhandle, "left_orientation"), "N"))
    stopifnot(identical(h5_read_vector(rhandle, "right_orientation"), "N"))

    Q <- as.matrix(reloadDelayedObject(rhandle, "left_seed", version=version, ...))
    Qty <- as.matrix(reloadDelayedObject(rhandle, "right_seed", version=version, ...))
    seed <- new("ResidualMatrixSeed", .matrix = .matrix, Q = Q, Qty = Qty, transposed = transposed)
    DelayedArray(seed)
}
