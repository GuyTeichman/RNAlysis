# Open a connection to a log file
logfile <- file("*PLACEHOLDERPATH*/logfile.log", open = "a")
# Redirect both output and messages to the file and console
sink(logfile, append = TRUE, split = TRUE)

require("DESeq2")

count_data <- read.table("tests/test_files/big_counted.csv", header=TRUE, sep= ",", row.names = 1)
design_matrix <- read.table("tests/test_files/test_design_matrix.csv", header=TRUE, sep= ",")
dds <- DESeqDataSetFromMatrix(count_data, design_matrix, ~ condition + replicate + factor1 + factor2)
dds_res <- DESeq(dds)
res <- results(dds_res, contrast=c('condition', 'cond2', 'cond1'), cooksCutoff=TRUE)
res_ordered <- res[order(res$padj),]
write.csv(as.data.frame(res_ordered),file="*PLACEHOLDERPATH*/DESeq2_condition_cond2_vs_cond1.csv")
dds_lrt <- DESeq(dds, test="LRT", reduced=~ condition + replicate + factor2)
lrt_res <- results(dds_lrt, cooksCutoff=TRUE)
lrt_res_ordered <- lrt_res[order(lrt_res$padj),]
write.csv(as.data.frame(lrt_res_ordered),file="*PLACEHOLDERPATH*/DESeq2_factor1_LRT.csv")
dds_lrt <- DESeq(dds, test="LRT", reduced=~ condition + replicate + factor1)
lrt_res <- results(dds_lrt, cooksCutoff=TRUE)
lrt_res_ordered <- lrt_res[order(lrt_res$padj),]
write.csv(as.data.frame(lrt_res_ordered),file="*PLACEHOLDERPATH*/DESeq2_factor2_LRT.csv")
# get session info
sessionInfo()
# Close the sink to stop redirecting output
sink()
