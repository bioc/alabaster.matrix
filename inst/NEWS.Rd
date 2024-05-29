\name{alabaster.matrix News}
\title{alabaster.matrix News}
\encoding{UTF-8}

\section{Version 1.6.0}{\itemize{
\item Support the SVT_SparseMatrix version 1 class definition in \code{saveObject}.
However, note that this was not implemented for the soft-deprecated \code{writeSparseMatrix}, which now errors if such objects are passed in.

\item Added a \code{extract_sparse_array()} method for the \code{WrapperArraySeed} class, for some future-proofing when the seeds eventually make the switch.
}}
