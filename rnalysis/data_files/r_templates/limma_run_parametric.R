require("limma")

design_matrix <- read.table("$DESIGN_MATRIX", header=TRUE, sep= ",")
$DEFINE_FACTORS
design <- model.matrix($FORMULA, design_matrix)
colnames(design)[1] <- "Intercept" # rename the intercept column from "(Intercept)" to "Intercept"
colnames(design) <- make.names(colnames(design)) # make coefficient names syntactically correct

count_data <- read.table("$COUNT_MATRIX", header=TRUE, sep= ",", row.names = 1)
voom_object <- voom(count_data, design, plot=FALSE, save.plot=TRUE)

$RANDOM_EFFECT_FIT
