\name{alabaster.matrix News}
\title{alabaster.matrix News}
\encoding{UTF-8}

\section{Version 1.14.0}{\itemize{
\item Add an explicit \code{saveObject()} method for HDF5Arrays, to avoid creating an unnecessary no-op file via the DelayedArray method.

\item Officially deprecate \code{stageObject()} and \code{loadObject()} in favor of \code{saveObject()} and \code{readObject()}.
This includes related functions like \code{writeSparseMatrix()}, \code{.createRawArraySeed()}, etc.
}}

\section{Version 1.8.0}{\itemize{
\item Deprecated the WrapperArraySeed in favor of direct inheritance from DelayedUnaryIsoOp.

\item Support deduplication of arrays and sparse matrices via the new \code{array.dedup.session=} argument in their corresponding \code{saveObject()} methods.
This also allows deduplication of external DelayedArray seeds if they have already been processed outside a DelayedArray in the same session.

\item Added a \code{DelayedArray.force.external=} argument to the \code{saveObject()} method for DelayedArrays.
This is a more convenient way of forcing seeds to be saved as external arrays for deduplication.

\item Support the new custom variable length string format for storing character arrays in HDF5.
This is controlled by the \code{array.character.vls=} option in the array's \code{saveObject()} method.
}}

\section{Version 1.6.0}{\itemize{
\item Support the SVT_SparseMatrix version 1 class definition in \code{saveObject()}.
However, note that this was not implemented for the soft-deprecated \code{writeSparseMatrix()}, which now errors if such objects are passed in.

\item Added a \code{extract_sparse_array()} method for the WrapperArraySeed class, for some future-proofing when the seeds eventually make the switch.

\item Bugfix for integer overflow when saving large sparse matrices in \code{saveObject()}.

\item Open all HDF5 files in read-only mode for \code{readObject()} dispatch functions, to avoid permission-related issues.

\item Added \code{altReloadDelayedObject()}, \code{altStoreDelayedObject()}, and their associated getters/setters,
to allow applications to override the delayed operation saving/reading process.

\item Added \code{registerReloadDelayedObjectFunction()} to allow extension developers to register reader functions for new classes.

\item Added a \code{ReloadedArray.reuse.files="relsymlink"} option in the \code{saveObject()} method for ReloadedArrays.
This creates relative symbolic links to the original array files, which is more robust to their movement provided the linked files are moved in the same manner.

\item Enable deduplication of identical seeds across multiple calls to \code{storeDelayedObject()} within a single \dQuote{session}.
This avoids making multiple copies of the same seed for different DelayedArray instances with the same seeds, e.g., in a SummarizedExperiment.

\item Added an \code{external.save.args=} option to \code{storeDelayedObject()} to avoid conflicts in the method arguments of \code{saveObject()}.
}}
