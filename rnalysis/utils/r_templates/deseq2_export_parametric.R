res <- results(dds_res, contrast=c$CONTRAST)
res_ordered <- res[order(res$padj),]
write.csv(as.data.frame(res_ordered),file="$OUTFILE_NAME")
