require("Rsubread")

fc <- featureCounts(arg = TRUE,
arg2 = NULL,
arg3 = "string",
arg4 = c("str1", "str2", "str 3"))

counts <- fc$counts
annotation <- fc$annotation
stats <- fc$stat

write.csv(as.data.frame(counts),file="C:/specified/output/dir/featureCounts_counts.csv")
write.csv(as.data.frame(annotation),file="C:/specified/output/dir/featureCounts_annotation.csv")
write.csv(as.data.frame(stats),file="C:/specified/output/dir/featureCounts_stats.csv")
