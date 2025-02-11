# Open a connection to a log file
logfile <- file("tests/test_files/fastq_pipeline_tests/single/outdir/03_featurecounts_single_end/logfile.log", open = "a")
# Redirect both output and messages to the file and console
sink(logfile, append = TRUE, split = TRUE)

require("Rsubread")

fc <- featureCounts(files = c("tests/test_files/fastq_pipeline_tests/single/outdir/02_bowtie2_align_single_end/reads_1_trimmed.sam", "tests/test_files/fastq_pipeline_tests/single/outdir/02_bowtie2_align_single_end/reads_2_trimmed.sam"),
annot.ext = "tests/test_files/kallisto_tests/transcripts.gtf",
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
nthreads = 1)

counts <- fc$counts
annotation <- fc$annotation
stats <- fc$stat

write.csv(as.data.frame(counts),file="tests/test_files/fastq_pipeline_tests/single/outdir/03_featurecounts_single_end/featureCounts_counts.csv")
write.csv(as.data.frame(annotation),file="tests/test_files/fastq_pipeline_tests/single/outdir/03_featurecounts_single_end/featureCounts_annotation.csv")
write.csv(as.data.frame(stats),file="tests/test_files/fastq_pipeline_tests/single/outdir/03_featurecounts_single_end/featureCounts_stats.csv")
# get session info
sessionInfo()
# Close the sink to stop redirecting output
sink()
