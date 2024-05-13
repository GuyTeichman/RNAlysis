dds_lrt <- DESeq(dds, test="LRT", reduced=$REDUCED)
lrt_res <- results(dds_lrt, cooksCutoff=$COOKS)
lrt_res_ordered <- lrt_res[order(lrt_res$padj),]
write.csv(as.data.frame(lrt_res_ordered),file="$OUTFILE_NAME")
