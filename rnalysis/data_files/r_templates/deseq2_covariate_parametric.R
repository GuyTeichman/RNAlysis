cov_res <- results(dds_res, name="$COVARIATE", cooksCutoff=$COOKS)
cov_res_ordered <- cov_res[order(res$padj),]
write.csv(as.data.frame(cov_res_ordered),file="$OUTFILE_NAME")
