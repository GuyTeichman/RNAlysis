# Open a connection to a log file
logfile <- file("*PLACEHOLDERPATH*/logfile.log", open = "a")
# Redirect both output and messages to the file and console
sink(logfile, append = TRUE, split = TRUE)

require("limma")

design_matrix <- read.table("tests/test_files/test_design_matrix.csv", header=TRUE, sep= ",")
design_matrix$condition <- factor(design_matrix$condition, levels=c("cond1", "cond2", "cond3"))
design_matrix$replicate <- factor(design_matrix$replicate, levels=c("rep1", "rep2", "rep3"))

design <- model.matrix(~ condition + replicate, design_matrix)
colnames(design)[1] <- "Intercept" # rename the intercept column from "(Intercept)" to "Intercept"
colnames(design) <- make.names(colnames(design)) # make coefficient names syntactically correct

count_data <- read.table("counted.csv", header=TRUE, sep= ",", row.names = 1)
voom_object <- voom(count_data, design, plot=FALSE, save.plot=TRUE)

fit <- lmFit(voom_object, design)

contrast <- makeContrasts("conditioncond3 - conditioncond2", levels = design)
contrast_fit <- contrasts.fit(fit, contrast)
contrast_bayes <- eBayes(contrast_fit)
res <- topTable(contrast_bayes, n=Inf)
res_ordered <- res[order(res$adj.P.Val),]
write.csv(as.data.frame(res_ordered),file="*PLACEHOLDERPATH*/LimmaVoom_condition_cond3_vs_cond2.csv")

contrast <- makeContrasts("replicaterep2", levels = design)
contrast_fit <- contrasts.fit(fit, contrast)
contrast_bayes <- eBayes(contrast_fit)
res <- topTable(contrast_bayes, n=Inf)
res_ordered <- res[order(res$adj.P.Val),]
write.csv(as.data.frame(res_ordered),file="*PLACEHOLDERPATH*/LimmaVoom_replicate_rep2_vs_rep1.csv")

contrast <- makeContrasts(" - conditioncond2", levels = design)
contrast_fit <- contrasts.fit(fit, contrast)
contrast_bayes <- eBayes(contrast_fit)
res <- topTable(contrast_bayes, n=Inf)
res_ordered <- res[order(res$adj.P.Val),]
write.csv(as.data.frame(res_ordered),file="*PLACEHOLDERPATH*/LimmaVoom_condition_cond1_vs_cond2.csv")
# get session info
sessionInfo()
# Close the sink to stop redirecting output
sink()
