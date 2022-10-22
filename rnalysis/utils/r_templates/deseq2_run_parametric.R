require("DESeq2")

count_data <- read.table("$COUNT_MATRIX", header=TRUE, sep= ",", row.names = 1)
design_matrix <- read.table("$DESIGN_MATRIX", header=TRUE, sep= ",")
dds <- DESeqDataSetFromMatrix(count_data, design_matrix, $FORMULA)
dds_res <- DESeq(dds)
