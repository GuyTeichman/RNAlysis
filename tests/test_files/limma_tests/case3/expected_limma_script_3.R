# Open a connection to a log file
logfile <- file("*PLACEHOLDERPATH*/logfile.log", open = "a")
# Redirect both output and messages to the file and console
sink(logfile, append = TRUE, split = TRUE)

require("limma")

design_matrix <- read.table("tests/test_files/test_design_matrix.csv", header=TRUE, sep= ",")
# name the design matrix rows after the samples (read as text, exactly as written in the file)
rownames(design_matrix) <- read.table("tests/test_files/test_design_matrix.csv", header=TRUE, sep= ",", colClasses="character", na.strings=character(0))[[1]]
design_matrix$condition <- factor(design_matrix$condition, levels=c("cond1", "cond2", "cond3"))
design_matrix$replicate <- factor(design_matrix$replicate, levels=c("rep1", "rep2", "rep3"))

design <- model.matrix(~ condition, design_matrix)
colnames(design)[1] <- "Intercept" # rename the intercept column from "(Intercept)" to "Intercept"
colnames(design) <- make.names(colnames(design)) # make coefficient names syntactically correct

count_data <- read.table("tests/test_files/big_counted.csv", header=TRUE, sep= ",", row.names = 1, check.names=FALSE)
# verify that the design matrix and the count matrix describe the same samples in the same order
if (!identical(rownames(design_matrix), colnames(count_data))) {
  stop("The sample names in the design matrix do not match the sample names in the count matrix")
}
voom_object <- voom(count_data, design, plot=FALSE, save.plot=TRUE)

cor <- duplicateCorrelation(voom_object, design, block=design_matrix$replicate)
if (cor$consensus.correlation > 0) { #only include random effect if the correlation is positive
  fit <- lmFit(voom_object, design, block=design_matrix$replicate, correlation=cor$consensus.correlation)
}else {   fit <- lmFit(voom_object, design)}

contrast <- makeContrasts("conditioncond2", levels = design)
contrast_fit <- contrasts.fit(fit, contrast)
contrast_bayes <- eBayes(contrast_fit)
res <- topTable(contrast_bayes, n=Inf)
res_ordered <- res[order(res$adj.P.Val),]
write.csv(as.data.frame(res_ordered),file="*PLACEHOLDERPATH*/LimmaVoom_condition_cond2_vs_cond1.csv")
# get session info
sessionInfo()
# Close the sink to stop redirecting output
sink()
