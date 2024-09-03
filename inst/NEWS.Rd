\name{alabaster.matrix News}
\title{alabaster.matrix News}
\encoding{UTF-8}

\section{Version 1.6.0}{\itemize{
\item Support the SVT_SparseMatrix version 1 class definition in \code{saveObject}.
However, note that this was not implemented for the soft-deprecated \code{writeSparseMatrix}, which now errors if such objects are passed in.

\item Added a \code{extract_sparse_array()} method for the \code{WrapperArraySeed} class, for some future-proofing when the seeds eventually make the switch.

\item Bugfix for integer overflow when saving large sparse matrices in \code{saveObject()}.

\item Open all HDF5 files in read-only mode for \code{readObject()} dispatch functions, to avoid permission-related issues.

\item Added \code{altReloadDelayedObject()}, \code{altStoreDelayedObject()}, and their associated getters/setters,
to allow applications to override the delayed operation saving/reading process.

\item Added \code{registerReloadDelayedObjectFunction()} to allow extension developers to register reader functions for new classes.

\item Added a \code{ReloadedArray.reuse.files="relsymlink"} option to create relative symbolic links to the original array files.
This is more robust to movement of the original files, provided the new files are moved in the same manner.

\item Enable deduplication of identical seeds across multiple calls to \code{storeDelayedObject()} within a single \dQuote{session}.
This avoids making multiple copies of the same seed for different DelayedArray instances with the same seeds, e.g., in a SummarizedExperiment.
}}
