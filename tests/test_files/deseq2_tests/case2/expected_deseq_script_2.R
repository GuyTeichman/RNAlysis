# Open a connection to a log file
logfile <- file("*PLACEHOLDERPATH*/logfile.log", open = "a")
# Redirect both output and messages to the file and console
sink(logfile, append = TRUE, split = TRUE)

require("DESeq2")

count_data <- read.table("counted.csv", header=TRUE, sep= ",", row.names = 1, check.names=FALSE)
design_matrix <- read.table("tests/test_files/test_design_matrix.csv", header=TRUE, sep= ",")
# name the design matrix rows after the samples (read as text, exactly as written in the file),
# so that R verifies that the design matrix and the count matrix describe the same samples in the same order
rownames(design_matrix) <- read.table("tests/test_files/test_design_matrix.csv", header=TRUE, sep= ",", colClasses="character", na.strings=character(0))[[1]]
if (!identical(rownames(design_matrix), colnames(count_data))) {
  stop("The sample names in the design matrix do not match the sample names in the count matrix")
}
dds <- DESeqDataSetFromMatrix(count_data, design_matrix, ~ condition + replicate)
dds_res <- DESeq(dds)
res <- results(dds_res, contrast=c("condition", "cond3", "cond2"), cooksCutoff=TRUE)
res_ordered <- res[order(res$padj),]
write.csv(as.data.frame(res_ordered),file="*PLACEHOLDERPATH*/DESeq2_condition_cond3_vs_cond2.csv")
res <- results(dds_res, contrast=c("replicate", "rep2", "rep1"), cooksCutoff=TRUE)
res_ordered <- res[order(res$padj),]
write.csv(as.data.frame(res_ordered),file="*PLACEHOLDERPATH*/DESeq2_replicate_rep2_vs_rep1.csv")
res <- results(dds_res, contrast=c("condition", "cond1", "cond2"), cooksCutoff=TRUE)
res_ordered <- res[order(res$padj),]
write.csv(as.data.frame(res_ordered),file="*PLACEHOLDERPATH*/DESeq2_condition_cond1_vs_cond2.csv")
# get session info
sessionInfo()
# Close the sink to stop redirecting output
sink()
