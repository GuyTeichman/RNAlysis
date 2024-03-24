require("limma")

design_matrix <- read.table("tests/test_files/test_design_matrix.csv", header=TRUE, sep= ",")
design_matrix$condition <- factor(design_matrix$condition, levels=c("cond1", "cond2", "cond3"))
design_matrix$replicate <- factor(design_matrix$replicate, levels=c("rep1", "rep2", "rep3"))

design <- model.matrix(~ condition + replicate, design_matrix)
colnames(design)[1] <- "Intercept" # rename the intercept column from "(Intercept)" to "Intercept"

count_data <- read.table("tests/test_files/big_counted.csv", header=TRUE, sep= ",", row.names = 1)
voom_object <- voom(count_data, design, plot=FALSE, save.plot=TRUE)

cor <- duplicateCorrelation(voom_object, design, block=replicate)
if (cor$consensus.correlation > 0) { #only include random effect if the correlation is positive
  fit <- lmFit(voom_object, design, block=replicate, correlation=cor$consensus.correlation)
}else {   fit <- lmFit(voom_object, design)}

contrast <- makeContrasts("conditioncond2", levels = design)
contrast_fit <- contrasts.fit(fit, contrast)
contrast_bayes <- eBayes(contrast_fit)
res <- topTable(contrast_bayes, n=Inf)
res_ordered <- res[order(res$adj.P.Val),]
write.csv(as.data.frame(res_ordered),file="*PLACEHOLDERPATH*/LimmaVoom_condition_cond2_vs_cond1.csv")
