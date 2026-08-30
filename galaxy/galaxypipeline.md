# Galaxy upstream pipeline

This document describes the upstream RNA-seq pipeline run in Galaxy (usegalaxy.org), from raw read download through Salmon quantification. The executable workflow is saved alongside this file as `sakaiensis-rnaseq-upstream.ga` and can be imported directly into Galaxy. This document is the human-readable companion: it explains what each step does, the settings that matter, and two problems encountered along the way.

The pipeline takes SRR accessions as input and produces one Salmon quantification file (`quant.sf`) per sample, which are the inputs to the downstream DESeq2 analysis in `scripts/deseq2_analysis.R`.

## Data

BioProject PRJNA1377810 (Novak and Gardner 2026). Nine samples, three conditions in biological triplicate: 702 complex medium, 702 with maltose, and 702 with PET plastic. Sequencing was 2x150 bp paired-end on an Illumina HiSeq 4000, averaging about 57 million reads per sample.

## Steps

The workflow has two branches that meet at Salmon: a reads branch (download, quality control, trimming) and a reference branch (genome download, extraction, header fix). Tool versions below are taken from the exported workflow file.

### Reads branch

**Faster Download and Extract Reads in FASTQ** (fasterq_dump v3.1.1+galaxy0). Downloads sequencing reads from the SRA and converts them from SRA format to FASTQ. Paired-end data produces two files per sample (forward and reverse), kept separate for later use by Salmon.

**Flatten collection** (used at several points). Galaxy groups outputs into nested collections. Flattening turns a nested collection into a flat list so that downstream tools and reports keep each sample labelled by its SRR name rather than collapsing samples together.

**FastQC** (fastqc v0.74+galaxy1). Reads through each FASTQ file and tallies quality statistics (per-base quality, GC content, duplication, adapter content). It changes nothing; it is purely descriptive. FastQC was run on a flattened collection so each of the 18 files kept its sample name.

**MultiQC** (multiqc v1.35+galaxy2). Aggregates the individual FastQC reports into one summary across all samples. MultiQC takes the FastQC raw data output, not the HTML report, as input.

**fastp** (fastp v1.3.6+galaxy0). Trims adapters and low-quality bases from the reads. Adapter detection found Nextera on read 1 and TruSeq on read 2; about 97 to 98% of reads were retained after trimming.

### Reference branch

**NCBI Datasets Genomes** (datasets_download_genome v18.33.1+galaxy0). Downloads the *P. sakaiensis* genome package from NCBI, including the coding sequence (CDS) FASTA and protein FASTA. The CDS FASTA is the reference used for quantification and, later, for gene annotation.

**Extract dataset**. Pulls the CDS FASTA out of the downloaded collection as a standalone dataset. This step is necessary because Salmon's transcriptome input requires a plain dataset, and Galaxy silently rejected the collection-wrapped form (see problem 1 below).

**Replace Text** (tp_replace_in_line v9.5+galaxy3). Replaces the pipe character in the FASTA headers. The find pattern is `\|` and the replacement is `_` (see problem 2 below). This produces the pipe-corrected CDS FASTA that the Salmon index is built from and that the annotation step in R parses.

### Quantification

**Salmon quant** (salmon v1.10.1+galaxy4). Builds an index from the pipe-corrected CDS FASTA and quantifies each sample against it. Salmon determines which gene each read came from and tallies the counts, without computing each read's exact position, which makes it fast and appropriate for a bacterium with no splicing. Library type was set to automatic and detected as IU (inward, unstranded). The index contained 5,523 targets after removing 49 duplicate sequences. Output is one `quant.sf` per sample with estimated read counts.

Mapping rates by group: complex about 71.7%, maltose about 73.9%, PET about 62.4%. All three PET samples mapped lower than every sugar-grown sample, consistent with the PET cultures being sampled later (about 27 hours) and nearer stationary phase, where a CDS-only index captures a smaller fraction of the transcriptome. Every sample still had ample assigned fragments for differential expression.

## Two problems encountered

**1. Collection-wrapped reference.** Salmon's transcriptome input needs a plain dataset, but the NCBI Datasets output is collection-wrapped. Galaxy allowed the collection to be selected in the form but silently rejected it on submission, with no error message and the form simply resetting. The fix was to use Extract dataset to pull the CDS FASTA out as a standalone dataset.

**2. Pipe characters in FASTA headers.** The CDS FASTA headers look like `>lcl|NZ_CP166677.1_cds_WP_054021747.1_1 [gene=dnaA] ...`. Salmon treats the pipe character as a name separator, so every gene was named `lcl` and indexing failed with an error about two references sharing the same name. The fix was the Replace Text step, replacing `\|` with `_`. This is why the CDS identifiers use underscores throughout the downstream analysis.

## Output

Nine `quant.sf` files, renamed by sample (`complex_1.sf` through `pet_3.sf`) and placed in `data/counts/` for the downstream DESeq2 analysis. The pipe-corrected CDS FASTA is saved as `data/reference/cds_reference.fasta` and is used both as the Salmon reference and as the source of gene annotation.