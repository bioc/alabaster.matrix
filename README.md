# Save array-like objects to file

|Environment|Status|
|---|---|
|[BioC-release](https://bioconductor.org/packages/release/bioc/html/alabaster.matrix.html)|[![Release OK](https://bioconductor.org/shields/build/release/bioc/alabaster.matrix.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/alabaster.matrix/)|
|[BioC-devel](https://bioconductor.org/packages/devel/bioc/html/alabaster.matrix.html)|[![Devel OK](https://bioconductor.org/shields/build/devel/bioc/alabaster.matrix.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/alabaster.matrix/)|

The **alabaster.matrix** package implements methods for saving and loading matrix- or array-like objects under the **alabaster** framework.
It provides a language-agnostic method for serializing data in arrays or abstractions thereof.
To get started, install the package and its dependencies from Bioconductor:

```r
# install.packages("BiocManager")
BiocManager::install("alabaster.matrix")
```

We can then save a variety of matrices and arrays to file.
For example, a sparse matrix can be saved to a HDF5 file in compressed sparse column format:

```r
library(Matrix)
y <- rsparsematrix(1000, 100, density=0.05)

# Saving it to a directory.
library(alabaster.matrix)
tmp <- tempfile()
saveObject(y, tmp)

# Reading it as a file-backed matrix.
roundtrip <- readObject(tmp)
roundtrip
## <1000 x 100> sparse ReloadedMatrix object of type "double":
##           [,1]   [,2]   [,3] ...  [,99] [,100]
##    [1,]      0      0      0   .      0      0
##    [2,]      0      0      0   .      0      0
##    [3,]      0      0      0   .      0      0
##    [4,]      0      0      0   .      0      0
##    [5,]      0      0      0   .      0      0
##     ...      .      .      .   .      .      .
##  [996,]      0      0      0   .      0      0
##  [997,]      0      0      0   .      0      0
##  [998,]      0      0      0   .      0      0
##  [999,]      0      0      0   .      0      0
## [1000,]      0      0      0   .      0      0

# Coerce this back into an in-memory sparse matrix:
inmemory <- as(roundtrip, "dgCMatrix")
```

We can also handle [`DelayedArray`](https://bioconductor.org/packages/DelayedArray) objects, possibly with preservation of delayed operations.
This uses the [**chihaya** specification](https://github.com/ArtifactDB/chihaya) to represent delayed operations inside a HDF5 file.

```r
library(DelayedArray)
y <- DelayedArray(rsparsematrix(1000, 100, 0.05))
y <- log1p(abs(y) / 1:100) # adding some delayed ops.

# Default method saves without preserving delayed operations.
tmp <- tempfile()
saveObject(y, tmp)
readObjectFile(tmp)$type
## [1] "compressed_sparse_matrix"

# But we can enable the delayed'ness explicitly, if so desired.
tmp2 <- tempfile()
saveObject(y, tmp2, delayedarray.preserve.ops=TRUE)
readObjectFile(tmp2)$type
## [1] "delayed_array"

roundtrip <- readObject(tmp2)
roundtrip
## <1000 x 100> sparse ReloadedMatrix object of type "double":
##           [,1]   [,2]   [,3] ...       [,99]      [,100]
##    [1,]      0      0      0   .           0           0
##    [2,]      0      0      0   .           0           0
##    [3,]      0      0      0   .           0           0
##    [4,]      0      0      0   .           0           0
##    [5,]      0      0      0   .           0           0
##     ...      .      .      .   .           .           .
##  [996,]      0      0      0   . 0.000000000 0.007368618
##  [997,]      0      0      0   . 0.000000000 0.000000000
##  [998,]      0      0      0   . 0.000000000 0.000000000
##  [999,]      0      0      0   . 0.000000000 0.000000000
## [1000,]      0      0      0   . 0.000000000 0.000000000
``` 
