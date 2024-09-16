# Open a connection to a log file
logfile <- file("*PLACEHOLDERPATH*/logfile.log", open = "a")
# Redirect both output and messages to the file and console
sink(logfile, append = TRUE, split = TRUE)

require("limma")

design_matrix <- read.table("tests/test_files/test_design_matrix_advanced.csv", header=TRUE, sep= ",")
design_matrix$condition <- factor(design_matrix$condition, levels=c("cond1", "cond2", "cond3"))
design_matrix$replicate <- factor(design_matrix$replicate, levels=c("rep1", "rep2", "rep3"))
design_matrix$covariate1 <- as.numeric(design_matrix$covariate1)
design_matrix$covariate2 <- as.numeric(design_matrix$covariate2)
design_matrix$factor1 <- factor(design_matrix$factor1, levels=c("A", "B", "C"))
design_matrix$factor2 <- factor(design_matrix$factor2, levels=c("X", "Y", "Z"))

design <- model.matrix(~ condition + replicate + covariate1 + covariate2 + factor1 + factor2 + factor1 + factor2, design_matrix)
colnames(design)[1] <- "Intercept" # rename the intercept column from "(Intercept)" to "Intercept"
colnames(design) <- make.names(colnames(design)) # make coefficient names syntactically correct

count_data <- read.table("tests/test_files/big_counted.csv", header=TRUE, sep= ",", row.names = 1)
voom_object <- voom(count_data, design, plot=FALSE, save.plot=TRUE)

fit <- lmFit(voom_object, design)

contrast <- makeContrasts("conditioncond2", levels = design)
contrast_fit <- contrasts.fit(fit, contrast)
contrast_bayes <- eBayes(contrast_fit)
res <- topTable(contrast_bayes, n=Inf)
res_ordered <- res[order(res$adj.P.Val),]
write.csv(as.data.frame(res_ordered),file="*PLACEHOLDERPATH*/LimmaVoom_condition_cond2_vs_cond1.csv")

bayes_fit <- eBayes(fit)
lrt_res <- topTable(bayes_fit, n=Inf, coef=c("factor1B", "factor1C"))
lrt_res_ordered <- lrt_res[order(lrt_res$adj.P.Val),]
write.csv(as.data.frame(lrt_res_ordered),file="*PLACEHOLDERPATH*/LimmaVoom_factor1_LRT.csv")

bayes_fit <- eBayes(fit)
lrt_res <- topTable(bayes_fit, n=Inf, coef=c("factor2Y", "factor2Z"))
lrt_res_ordered <- lrt_res[order(lrt_res$adj.P.Val),]
write.csv(as.data.frame(lrt_res_ordered),file="*PLACEHOLDERPATH*/LimmaVoom_factor2_LRT.csv")
# get session info
sessionInfo()
# Close the sink to stop redirecting output
sink()
