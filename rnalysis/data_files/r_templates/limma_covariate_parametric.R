
bayes_fit <- eBayes(fit)
cov_res <- topTable(bayes_fit, n=Inf, coef=$COEF)
cov_res_ordered <- cov_res[order(res$adj.P.Val),]
write.csv(as.data.frame(cov_res_ordered),file="$OUTFILE_NAME")
