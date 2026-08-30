# plastic-enzyme-rnaseq

RNA-seq differential expression analysis of *Piscinibacter sakaiensis* 201-F6 grown on PET plastic versus simple sugars.

## Question

When *P. sakaiensis* is grown on PET plastic rather than a simple sugar, which genes does it turn on, and are the known plastic-degradation enzymes (PETase and MHETase) among them?

## Data

Raw reads from BioProject PRJNA1377810 (Novak and Gardner 2026): nine samples across three conditions (702 complex medium, 702 with maltose, 702 with PET plastic) in biological triplicate. The source paper published the data sets only and did not perform differential expression analysis.

## Key finding

PETase is strongly induced on PET (roughly 29-fold), and MHETase is also significantly upregulated (roughly two-fold). Downstream ethylene glycol and aromatic degradation genes show smaller and more mixed changes. See `report/report.md` for the full write-up and limitations.

## Pipeline

Upstream steps (read download, quality control, trimming, and Salmon quantification) were run in Galaxy. Downstream differential expression was run in R with DESeq2. Salmon quantification files (`data/counts/`) and the coding sequence reference (`data/reference/`) are included so the R analysis runs from a clean clone. Raw reads and intermediate Galaxy files are excluded.

## Repository structure

```
data/
  counts/       salmon quantification files (.sf), one per sample
  metadata/     sample table (condition, replicate, accessions)
  reference/    coding sequence fasta used for quantification and annotation
results/
  de/           differential expression tables and gene annotation
  figures/      pca plot and pathway gene figure
scripts/
  deseq2_analysis.R    full differential expression pipeline
  pathway_figure.R     pathway gene bar chart
report/
  report.md     results, interpretation, and limitations
```

## How to run

1. Open the repository in R with the working directory set to the repository root.
2. Run `scripts/deseq2_analysis.R` to reproduce the differential expression tables and PCA plot.
3. Run `scripts/pathway_figure.R` to reproduce the pathway gene figure.

## References

Novak, J.K., and Gardner, J.G. 2026. Transcriptome data sets of *Piscinibacter sakaiensis* grown in 702 complex medium, maltose, and PET plastic. Microbiology Resource Announcements. 15(4).

Yoshida, S., et al. 2016. A bacterium that degrades and assimilates poly(ethylene terephthalate). Science. 351: 1196-1199.