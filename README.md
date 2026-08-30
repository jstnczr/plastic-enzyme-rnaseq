# plastic-enzyme-rnaseq

RNA-seq differential expression analysis of *Piscinibacter sakaiensis* 201-F6 grown on PET plastic versus simple sugars. Built in Galaxy (upstream) and R (downstream).

---

## Question

When *P. sakaiensis* is grown on PET plastic rather than a simple sugar, which genes does it turn on, and are the known plastic-degradation enzymes (PETase and MHETase) among them?

## Key finding

PETase is strongly induced on PET (roughly 29-fold), and MHETase is also significantly upregulated (roughly two-fold). Downstream ethylene glycol and aromatic degradation genes show smaller and more mixed changes. See [`report/report.md`](report/report.md) for the full write-up and limitations.

---

## Data

Raw reads from BioProject PRJNA1377810 (Novak and Gardner 2026): nine samples across three conditions (702 complex medium, 702 with maltose, 702 with PET plastic) in biological triplicate. The source paper published the data sets only and did not perform differential expression analysis.

## Pipeline

Upstream steps (read download, quality control, trimming, and Salmon quantification) were run in Galaxy and are documented in [`galaxy/`](galaxy/), including the importable workflow file. Downstream differential expression was run in R with DESeq2. Salmon quantification files (`data/counts/`) and the coding sequence reference (`data/reference/`) are included so the R analysis runs from a clean clone. Raw reads and intermediate Galaxy files are excluded.

---

## Repository Structure

```
plastic-enzyme-rnaseq/
├── data/
│   ├── counts/       salmon quantification files (.sf), one per sample
│   ├── metadata/     sample table (condition, replicate, accessions)
│   └── reference/    coding sequence fasta used for quantification and annotation
├── results/
│   ├── de/           differential expression tables and gene annotation
│   └── figures/      pca plot and pathway gene figure
├── scripts/
│   ├── deseq2_analysis.R    full differential expression pipeline
│   └── pathway_figure.R     pathway gene bar chart
├── galaxy/
│   ├── sakaiensis-rnaseq-upstream.ga   importable galaxy workflow
│   └── pipeline.md                     human-readable pipeline description
└── report/
    └── report.md     results, interpretation, and limitations
```

---

## How to Reproduce

Requirements: R with DESeq2, tximport, apeglm, and Biostrings (from Bioconductor), plus ggplot2 (from CRAN). The upstream Galaxy workflow is provided separately in `galaxy/` and does not need to be rerun, since its outputs (the `.sf` files and reference) are included.

```r
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "tximport", "apeglm", "Biostrings"))
install.packages("ggplot2")
```

With the working directory set to the repository root:

```r
# reproduce the differential expression tables and pca plot
source("scripts/deseq2_analysis.R")

# reproduce the pathway gene figure
source("scripts/pathway_figure.R")
```

---

## Tools

Galaxy (fasterq-dump, FastQC, MultiQC, fastp, NCBI Datasets, Salmon), R, DESeq2, tximport, apeglm, Biostrings, ggplot2

---

## Background

This project is a follow-up to an earlier comparative and structural analysis of plastic-degrading enzymes, extending that work from sequence and structure to gene expression: rather than asking what the PET-degradation enzymes look like, it asks when the organism actually turns them on. *P. sakaiensis* is the bacterium in which PETase and MHETase were first characterized (Yoshida et al. 2016), making it a natural system for asking how the plastic-degradation pathway is regulated.

---

## References

Novak, J.K., and Gardner, J.G. 2026. Transcriptome data sets of *Piscinibacter sakaiensis* grown in 702 complex medium, maltose, and PET plastic. Microbiology Resource Announcements. 15: e01466-25. https://doi.org/10.1128/mra.01466-25

Yoshida, S., et al. 2016. A bacterium that degrades and assimilates poly(ethylene terephthalate). Science. 351: 1196-1199. https://doi.org/10.1126/science.aad6359