
contrast <- makeContrasts($CONTRAST, levels = design)
contrast_fit <- contrasts.fit(fit, contrast)
contrast_bayes <- eBayes(contrast_fit)
res <- topTable(contrast_bayes, n=Inf)
res_ordered <- res[order(res$adj.P.Val),]
write.csv(as.data.frame(res_ordered),file="$OUTFILE_NAME")
