require("limma")

design_matrix <- read.table("$DESIGN_MATRIX", header=TRUE, sep= ",")
$DEFINE_FACTORS
design <- model.matrix($FORMULA)
colnames(design)[1] <- "Intercept"

count_data <- read.table("$COUNT_MATRIX", header=TRUE, sep= ",", row.names = 1)
voom_object <- voom(count_data, design, plot=FALSE, save.plot=TRUE)

$RANDOM_EFFECT_FIT
