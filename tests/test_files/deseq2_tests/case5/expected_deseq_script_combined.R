require("DESeq2")

count_data <- read.table("tests/test_files/big_counted.csv", header=TRUE, sep= ",", row.names = 1)
design_matrix <- read.table("tests/test_files/test_design_matrix.csv", header=TRUE, sep= ",")
dds <- DESeqDataSetFromMatrix(count_data, design_matrix, ~ condition + replicate + covariate1 + covariate2 + factor1 + factor2)
dds_res <- DESeq(dds)
res <- results(dds_res, contrast=c('condition', 'cond2', 'cond1'))
res_ordered <- res[order(res$padj),]
write.csv(as.data.frame(res_ordered),file="*PLACEHOLDERPATH*/DESeq2_condition_cond2_vs_cond1.csv")
cov_res <- results(dds_res, name="covariate1")
cov_res_ordered <- cov_res[order(res$padj),]
write.csv(as.data.frame(cov_res_ordered),file="*PLACEHOLDERPATH*/DESeq2_covariate1_covariate.csv")
cov_res <- results(dds_res, name="covariate2")
cov_res_ordered <- cov_res[order(res$padj),]
write.csv(as.data.frame(cov_res_ordered),file="*PLACEHOLDERPATH*/DESeq2_covariate2_covariate.csv")
dds_lrt <- DESeq(dds, test="LRT", reduced=~ condition + replicate + covariate1 + covariate2 + factor2)
resultsNames(dds_lrt)
lrt_res <- results(dds_lrt)
lrt_res_ordered <- lrt_res[order(lrt_res$padj),]
write.csv(as.data.frame(lrt_res_ordered),file="*PLACEHOLDERPATH*/DESeq2_factor1_LRT.csv")
dds_lrt <- DESeq(dds, test="LRT", reduced=~ condition + replicate + covariate1 + covariate2 + factor1)
resultsNames(dds_lrt)
lrt_res <- results(dds_lrt)
lrt_res_ordered <- lrt_res[order(lrt_res$padj),]
write.csv(as.data.frame(lrt_res_ordered),file="*PLACEHOLDERPATH*/DESeq2_factor2_LRT.csv")
