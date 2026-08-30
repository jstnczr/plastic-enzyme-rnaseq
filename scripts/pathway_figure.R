# figure: pet degradation pathway gene expression changes
# horizontal bar chart of named pathway enzymes, pet vs maltose
library(ggplot2)

# read the annotated results
de_annotated <- read.csv("results/de/pet_vs_maltose_annotated.csv")

# pick the named pathway enzymes we can tie to source papers
# using protein description text to select them
pathway_terms <- c("poly\\(ethylene terephthalate\\) hydrolase",
                   "mono\\(2-hydroxyethyl\\) terephthalate hydrolase",
                   "ethylene glycol dehydrogenase",
                   "protocatechuate 3,4-dioxygenase",
                   "glycolate oxidase",
                   "glyoxylate carboligase")

# build the regex and pull matching rows
pattern <- paste(pathway_terms, collapse = "|")
fig_data <- de_annotated[grepl(pattern, de_annotated$protein, ignore.case = TRUE), ]

# keep only significant genes (padj < 0.05)
fig_data <- fig_data[!is.na(fig_data$padj) & fig_data$padj < 0.05, ]

# make short readable labels instead of long protein names
fig_data$label <- fig_data$protein
fig_data$label <- gsub("poly\\(ethylene terephthalate\\) hydrolase", "PETase", fig_data$label)
fig_data$label <- gsub("mono\\(2-hydroxyethyl\\) terephthalate hydrolase", "MHETase", fig_data$label)
fig_data$label <- gsub(" subunit", "", fig_data$label)

# direction of change for colouring
fig_data$direction <- ifelse(fig_data$log2FoldChange > 0, "Up on PET", "Down on PET")

# order bars by fold change
fig_data <- fig_data[order(fig_data$log2FoldChange), ]
fig_data$label <- factor(fig_data$label, levels = fig_data$label)

# build the plot
p <- ggplot(fig_data, aes(x = log2FoldChange, y = label, fill = direction)) +
  geom_col() +
  scale_fill_manual(values = c("Up on PET" = "#2c7fb8", "Down on PET" = "#d95f0e")) +
  labs(title = "PET-degradation pathway genes respond to growth on PET",
       subtitle = "Piscinibacter sakaiensis, PET vs maltose (significant genes, padj < 0.05)",
       x = "log2 fold change (PET vs maltose)",
       y = NULL,
       fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# save it
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("results/figures/pathway_genes.png", p, width = 9, height = 5, dpi = 150)
print("saved results/figures/pathway_genes.png")

# also print the data so you can see exactly what's plotted
print(fig_data[, c("label", "log2FoldChange", "padj")])