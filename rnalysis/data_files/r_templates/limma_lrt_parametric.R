
bayes_fit <- eBayes(fit)
lrt_res <- topTable(bayes_fit, n=Inf, coef=$COEFS)
lrt_res_ordered <- lrt_res[order(res$adj.P.Val),]
write.csv(as.data.frame(lrt_res_ordered),file="$OUTFILE_NAME")
