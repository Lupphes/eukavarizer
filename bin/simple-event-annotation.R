#!/usr/bin/env Rscript

# Create a local R library folder if it doesn’t exist
dir.create("Rlibs", showWarnings = FALSE)
.libPaths("Rlibs")

# Install packages to the local lib
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
# BiocManager::install(c("VariantAnnotation", "StructuralVariantAnnotation"), ask = FALSE, update = FALSE)
# install.packages("stringr", repos = "https://cloud.r-project.org")

# Load libraries
library(VariantAnnotation)
library(StructuralVariantAnnotation)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_vcf <- args[2]
output_bed <- args[3]

vcf <- readVcf(input_file, "hg19")

info(header(vcf)) <- unique(as(rbind(
  as.data.frame(info(header(vcf))),
  data.frame(
    row.names = c("SIMPLE_TYPE", "SVLEN"),
    Number = c("1", "1"),
    Type = c("String", "Integer"),
    Description = c(
      "Simple event type annotation based purely on breakend position and orientation.",
      "Structural variant length"
    )
  )
), "DataFrame"))

simpleEventType <- function(gr) {
  pgr = partner(gr)
  return(ifelse(seqnames(gr) != seqnames(pgr), "CTX",
    ifelse(strand(gr) == strand(pgr), "INV",
      ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS",
        ifelse(xor(start(gr) < start(pgr), strand(gr) == "-"), "DEL", "DUP")))))
}

gr <- breakpointRanges(vcf)
svtype <- simpleEventType(gr)

info(vcf)$SIMPLE_TYPE <- NA_character_
info(vcf)$SVLEN <- NA_integer_

info(vcf)[gr$sourceId, "SIMPLE_TYPE"] <- svtype
info(vcf)[gr$sourceId, "SVLEN"] <- gr$svLen

writeVcf(vcf, output_vcf)

# Filter for PASS and simple SV types
gr <- gr[gr$FILTER == "PASS" & partner(gr)$FILTER == "PASS"]
simplegr <- gr[simpleEventType(gr) %in% c("INS", "INV", "DEL", "DUP")]

simplebed <- data.frame(
  chrom = seqnames(simplegr),
  start = as.integer((start(simplegr) + end(simplegr)) / 2),
  end   = as.integer((start(partner(simplegr)) + end(partner(simplegr))) / 2),
  name  = simpleEventType(simplegr),
  score = simplegr$QUAL,
  strand = "."
)

simplebed <- simplebed[simplebed$start < simplebed$end,]
write.table(simplebed, output_bed, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
