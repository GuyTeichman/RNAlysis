require("limma")

design_matrix <- read.table($DESIGN_MATRIX, header=TRUE, sep= ",")
# name the design matrix rows after the samples (read as text, exactly as written in the file)
rownames(design_matrix) <- read.table($DESIGN_MATRIX, header=TRUE, sep= ",", colClasses="character", na.strings=character(0))[[1]]
$DEFINE_FACTORS
design <- model.matrix($FORMULA, design_matrix)
colnames(design)[1] <- "Intercept" # rename the intercept column from "(Intercept)" to "Intercept"
colnames(design) <- make.names(colnames(design)) # make coefficient names syntactically correct

count_data <- read.table($COUNT_MATRIX, header=TRUE, sep= ",", row.names = 1, check.names=FALSE)
# verify that the design matrix and the count matrix describe the same samples in the same order
if (!identical(rownames(design_matrix), colnames(count_data))) {
  stop("The sample names in the design matrix do not match the sample names in the count matrix")
}
voom_object <- $VOOMFUNC(count_data, design, plot=FALSE, save.plot=TRUE)

$RANDOM_EFFECT_FIT
