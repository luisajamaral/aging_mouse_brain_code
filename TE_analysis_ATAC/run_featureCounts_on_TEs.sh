#!/bin/bash

# TE quantification using featureCounts
# Input: BAM files from split cell types, SAF annotation of TE subfamilies
# Output: Count matrix of TE subfamilies across samples

# Define SAF file and output
ANNOTATION="out.saf"
OUTPUT="allTE.feature.counts.txt"

# List of BAMs
BAMS=(*.bam)

# Run featureCounts
featureCounts \
  -a "$ANNOTATION" \
  -o "$OUTPUT" \
  -F SAF \
  -M -O \
  -T 16 \
  -p --countReadPairs \
  "${BAMS[@]}"