require("limma")

design_matrix <- read.table("tests/test_files/test_design_matrix.csv", header=TRUE, sep= ",")
condition <- factor(design_matrix$condition, levels=c("cond1", "cond2", "cond3"))
replicate <- factor(design_matrix$replicate, levels=c("rep1", "rep2", "rep3"))

design <- model.matrix(~ condition + replicate)
colnames(design)[1] <- "Intercept"

count_data <- read.table("tests/test_files/big_counted.csv", header=TRUE, sep= ",", row.names = 1)
voom_object <- voom(count_data, design, plot=FALSE, save.plot=TRUE)
fit <- lmFit(voom_object, design)


contrast <- makeContrasts("conditioncond2", levels = design)
contrast_fit <- contrasts.fit(fit, contrast)
contrast_bayes <- eBayes(contrast_fit)
res <- topTable(contrast_bayes, n=Inf)
res_ordered <- res[order(res$adj.P.Val),]
write.csv(as.data.frame(res_ordered),file="*PLACEHOLDERPATH*/LimmaVoom_condition_cond2_vs_cond1.csv")
