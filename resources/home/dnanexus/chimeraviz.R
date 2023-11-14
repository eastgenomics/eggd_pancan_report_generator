# this script is to make png of each plot per fusion.
# this is because chimeraviz takes up a lot of space so this script
# is designed to plot ONE fusion PER sample. It expects certain
# argument inputs in a specific order (in the working version this will have 
# parameter inputs set)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
suppressMessages(library(chimeraviz))
# args 1 = ensb sqlite
# args 2 = sample BAM
# args 3 = fusion_inspector table output
# args 4 = fusion_number


# get ensb annotation sqlite input from cmd line
edb <- ensembldb::EnsDb(args[1])
# get bam file input from cmd line
fusion_reads <- args[2]
# get fusion number input from cmd line
FI_fusions <- import_starfusion(args[3], "hg38")
fusion <- get_fusion_by_id(FI_fusions, args[4]) 
fusion_df <- read.csv(args[3], sep = "\t")
filename <- fusion_df[args[4],1]

## plot 1 - overall gene fusion
png(paste(filename, "gene_plot.png", sep = "_"))
plot_fusion_transcript(
    fusion,
    edb)
dev.off()
## plot 2 - per transcript with coverage
png(paste(filename, "transcript_plot.png", sep = "_"), width = 1000)
plot_fusion(
    fusion = fusion,
    bamfile = fusion_reads,
    edb = edb,
    non_ucsc = TRUE)
dev.off()

