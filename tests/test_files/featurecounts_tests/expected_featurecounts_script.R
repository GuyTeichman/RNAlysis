# Open a connection to a log file
logfile <- file("tests/test_files/featurecounts_tests/outdir/logfile.log", open = "a")
# Redirect both output and messages to the file and console
sink(logfile, append = TRUE, split = TRUE)

require("Rsubread")

fc <- featureCounts(arg = TRUE,
arg2 = NULL,
arg3 = "string",
arg4 = c("str1", "str2", "str 3"))

counts <- fc$counts
annotation <- fc$annotation
stats <- fc$stat

write.csv(as.data.frame(counts),file="tests/test_files/featurecounts_tests/outdir/featureCounts_counts.csv")
write.csv(as.data.frame(annotation),file="tests/test_files/featurecounts_tests/outdir/featureCounts_annotation.csv")
write.csv(as.data.frame(stats),file="tests/test_files/featurecounts_tests/outdir/featureCounts_stats.csv")
# get session info
sessionInfo()
# Close the sink to stop redirecting output
sink()
