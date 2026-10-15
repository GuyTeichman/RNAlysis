require("DESeq2")

count_data <- read.table($COUNT_MATRIX, header=TRUE, sep= ",", row.names = 1, check.names=FALSE)
design_matrix <- read.table($DESIGN_MATRIX, header=TRUE, sep= ",")
# name the design matrix rows after the samples (read as text, exactly as written in the file),
# so that R verifies that the design matrix and the count matrix describe the same samples in the same order
rownames(design_matrix) <- read.table($DESIGN_MATRIX, header=TRUE, sep= ",", colClasses="character", na.strings=character(0))[[1]]
if (!identical(rownames(design_matrix), colnames(count_data))) {
  stop("The sample names in the design matrix do not match the sample names in the count matrix")
}
dds <- DESeqDataSetFromMatrix(count_data, design_matrix, $FORMULA)
$NORMFACTORS
dds_res <- DESeq(dds)
