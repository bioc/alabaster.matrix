.onLoad <- function(libname, pkgname) {
    registerReadObjectFunction("dense_array", readArray)
    registerReadObjectFunction("compressed_sparse_matrix", readSparseMatrix)
    registerReadObjectFunction("delayed_array", readDelayedArray)
}

.onUnload <- function(libname, pkgname) {
    registerReadObjectFunction("dense_array", NULL)
    registerReadObjectFunction("compressed_sparse_matrix", NULL)
    registerReadObjectFunction("delayed_array", NULL)
}
