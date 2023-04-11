require("Rsubread")

fc <- featureCounts(files = c("tests/test_files/featurecounts_tests/paired/test-minimum.bam"),
annot.ext = "tests/test_files/featurecounts_tests/test-minimum.gtf",
isGTFAnnotationFile = TRUE,
GTF.featureType = "exon",
GTF.attrType = "gene_id",
allowMultiOverlap = FALSE,
countMultiMappingReads = FALSE,
fraction = FALSE,
isLongRead = FALSE,
minMQS = 0,
primaryOnly = TRUE,
strandSpecific = 0,
nthreads = 1,
isPairedEnd = TRUE,
requireBothEndsMapped = TRUE,
countChimericFragments = FALSE,
minFragLength = 50,
maxFragLength = 600,
countReadPairs = TRUE)

counts <- fc$counts
annotation <- fc$annotation
stats <- fc$stat

write.csv(as.data.frame(counts),file="tests/test_files/featurecounts_tests/outdir/featureCounts_counts.csv")
write.csv(as.data.frame(annotation),file="tests/test_files/featurecounts_tests/outdir/featureCounts_annotation.csv")
write.csv(as.data.frame(stats),file="tests/test_files/featurecounts_tests/outdir/featureCounts_stats.csv")
