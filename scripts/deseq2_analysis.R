# deseq2 differential expression analysis
# pet vs maltose in piscinibacter sakaiensis
# data: novak & gardner 2026, bioproject prjna1377810

# load packages
library(DESeq2)
library(tximport)
library(apeglm)

# read the sample metadata
samples <- read.csv("data/metadata/samples.csv")
print("loaded sample metadata")
print(samples)

# build paths to each .sf file from the sample names
files <- file.path("data", "counts", paste0(samples$sample_name, ".sf"))

# name the vector so columns get labelled by sample
names(files) <- samples$sample_name

# check every file actually exists before importing
print("checking that all .sf files exist...")
print(file.exists(files))

# import salmon quant files at transcript (cds) level, skipping inferential replicates
# txOut = TRUE because this is a bacterium, no splicing, one cds = one gene
txi <- tximport(files, type = "salmon", txOut = TRUE, dropInfReps = TRUE)
print("imported salmon quant files")

# turn condition into a factor and set maltose as the reference level
# maltose is the carbon-matched baseline, so pet_vs_maltose shrinks directly
samples$condition <- factor(samples$condition, levels = c("maltose", "complex", "pet"))

# build the deseq2 dataset from the tximport object
dds <- DESeqDataSetFromTximport(txi,
                                colData = samples,
                                design = ~ condition)
print("built DESeqDataSet")

# keep only genes with at least 10 total reads across all samples
keep <- rowSums(counts(dds)) >= 10
print("genes passing the count filter:")
print(table(keep))

# apply the filter
dds <- dds[keep, ]

# run the full deseq2 analysis (size factors, dispersion, glm fit)
dds <- DESeq(dds)
print("ran DESeq()")

# show the coefficient names that were estimated
print("coefficients estimated:")
print(resultsNames(dds))

# variance-stabilizing transform for the pca plot
vsd <- vst(dds, blind = TRUE)

# save the pca plot to the figures folder
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
png("results/figures/pca_plot.png", width = 1200, height = 700, res = 150)
plotPCA(vsd, intgroup = "condition")
dev.off()
print("saved pca plot to results/figures/pca_plot.png")

# extract pet vs maltose results with a stricter fold-change and significance threshold
# lfcThreshold = 1 tests for at least a 2-fold change, alpha = 0.05 sets 5% fdr
res_pet <- results(dds,
                   contrast = c("condition", "pet", "maltose"),
                   lfcThreshold = 1,
                   alpha = 0.05)
print("pet vs maltose summary at lfcThreshold = 1, alpha = 0.05:")
summary(res_pet)

# shrink the log2 fold changes for the pet vs maltose coefficient
# apeglm needs a named coefficient from resultsNames(dds)
res_pet_shrunk <- lfcShrink(dds,
                            coef = "condition_pet_vs_maltose",
                            type = "apeglm")
print("applied apeglm shrinkage")

# order by adjusted p-value
res_ordered <- res_pet_shrunk[order(res_pet_shrunk$padj), ]

# make sure the output directory exists
dir.create("results/de", recursive = TRUE, showWarnings = FALSE)

# convert results to a data frame and keep the cds id as a column
res_df <- as.data.frame(res_ordered)
res_df$cds_id <- rownames(res_df)

# reorder so cds_id is the first column
res_df <- res_df[, c("cds_id", setdiff(names(res_df), "cds_id"))]

# write the full shrunk results
write.csv(res_df, "results/de/pet_vs_maltose_deseq2.csv", row.names = FALSE)
print(paste("wrote", nrow(res_df), "rows to results/de/pet_vs_maltose_deseq2.csv"))

# --- annotation from the cds fasta headers ---
library(Biostrings)

# read the cds fasta
fasta <- readDNAStringSet("data/reference/cds_reference.fasta")

# the header is stored in the names of the object
headers <- names(fasta)
print(paste("read", length(headers), "headers"))

# cds_id is everything before the first space
cds_id <- sub(" .*", "", headers)

# pull each bracketed field with regex, returning NA if the field is absent
gene <- sub(".*\\[gene=([^]]+)\\].*", "\\1", headers)
gene[!grepl("\\[gene=", headers)] <- NA

locus_tag <- sub(".*\\[locus_tag=([^]]+)\\].*", "\\1", headers)
locus_tag[!grepl("\\[locus_tag=", headers)] <- NA

protein <- sub(".*\\[protein=([^]]+)\\].*", "\\1", headers)
protein[!grepl("\\[protein=", headers)] <- NA

protein_id <- sub(".*\\[protein_id=([^]]+)\\].*", "\\1", headers)
protein_id[!grepl("\\[protein_id=", headers)] <- NA

# assemble into a data frame
annotation <- data.frame(cds_id, gene, locus_tag, protein, protein_id,
                         stringsAsFactors = FALSE)

# write the annotation table
write.csv(annotation, "results/de/annotation.csv", row.names = FALSE)
print(paste("wrote", nrow(annotation), "annotation rows"))

# --- join annotation onto the de results ---

# read the de results back in
de <- read.csv("results/de/pet_vs_maltose_deseq2.csv")

# left join annotation onto the de results by cds_id
de_annotated <- merge(de, annotation, by = "cds_id", all.x = TRUE)

# re-sort by adjusted p-value since merge reorders rows
de_annotated <- de_annotated[order(de_annotated$padj), ]

# write the annotated results
write.csv(de_annotated, "results/de/pet_vs_maltose_annotated.csv", row.names = FALSE)
print(paste("wrote", nrow(de_annotated), "annotated rows"))

# --- targeted pet degradation pathway gene table ---

# search terms drawn from the pet degradation pathway (yoshida 2016)
# each is a known step in the pet degradation pathway
search_terms <- c("terephthalate",
                  "poly(ethylene terephthalate)",
                  "MHET",
                  "protocatechuate",
                  "ethylene glycol",
                  "glycolate",
                  "glyoxylate",
                  "aldehyde dehydrogenase")

# build a single regex that matches any of the terms
pattern <- paste(search_terms, collapse = "|")

# find annotation rows whose protein description matches any term
hits <- de_annotated[grepl(pattern, de_annotated$protein, ignore.case = TRUE), ]

# show the useful columns, sorted by fold change
hits_sorted <- hits[order(-hits$log2FoldChange), c("gene", "protein", "log2FoldChange", "padj")]

# write the targeted pathway gene table
write.csv(hits_sorted, "results/de/pet_pathway_genes.csv", row.names = TRUE)
print(paste("wrote", nrow(hits_sorted), "pathway gene rows"))

print("analysis complete")
