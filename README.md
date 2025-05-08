[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/mdusFg9r)

# microDNA Detection
Detection and validation of candidate microDNA circles using soft-clipped reads from aligned sequencing data.

## Overview
This pipeline identifies potential microDNA elements by detecting soft-clipped reads from a BAM file, clustering junction candidates, aligning clipped sequences, and scoring circle confidence based on microhomology and genomic distance. Final candidates are filtered based on read support and confidence thresholds.

---

## Usage

### Detect Soft-Clipped Reads
```bash
usage: detecting.py --bam BAM_FILE [--min_clip MIN_CLIP] [--out OUT_FILE]

optional arguments:
  --bam BAM_FILE     Input BAM file (indexed)
  --min_clip         Minimum clip length (default: 10)
  --out              Output TSV file for clipped reads (default: soft_clipped_reads.tsv)

---

## Validate Circles 

```bash

usage: validate_circles.py --circles CIRCLES_FILE --reads READS_FILE --start_positions START [START ...]

optional arguments:
  --circles          TSV file with candidate microDNA circles
  --reads            TSV file with clipped reads from detection step
  --start_positions  One or more start positions to validate (space-separated)

## Example: Validating Circle Candidates
```bash

$ python validate_circles.py \
    --circles microdna_circles.tsv \
    --reads clips_12bp.tsv \
    --start_positions 121485177 121485185 17912842

