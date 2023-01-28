require("Rsubread")

fc <- featureCounts($KWARGS)

counts <- fc$counts
annotation <- fc$annotation
stats <- fc$stat

write.csv(as.data.frame(counts),file="$OUTPUT_DIR/featureCounts_counts.csv")
write.csv(as.data.frame(annotation),file="$OUTPUT_DIR/featureCounts_annotation.csv")
write.csv(as.data.frame(stats),file="$OUTPUT_DIR/featureCounts_stats.csv")
