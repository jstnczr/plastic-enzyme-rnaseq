# Transcriptional response of *Piscinibacter sakaiensis* to growth on PET plastic

## Background

*Piscinibacter sakaiensis* 201-F6 (formerly *Ideonella sakaiensis*) is a gram-negative bacterium able to degrade and assimilate poly(ethylene terephthalate) (PET) plastic (Yoshida et al. 2016). Two enzymes carry out the initial breakdown of PET: PET hydrolase (PETase), which converts PET to mono(2-hydroxyethyl) terephthalate (MHET), and MHET hydrolase (MHETase), which converts MHET to terephthalate and ethylene glycol (Yoshida et al. 2016). Novak and Gardner (2026) published raw transcriptome data sets of *P. sakaiensis* grown on three carbon sources but did not analyze differential expression. This project uses those data to ask which genes respond when the organism is grown on PET rather than a simple sugar, and whether the known PET-degradation genes are among them.

## Methods

Raw reads were obtained from BioProject PRJNA1377810 (Novak and Gardner 2026): nine samples across three conditions (702 complex medium, 702 with maltose, 702 with PET plastic) in biological triplicate. Reads were quality-checked, trimmed, and quantified against the *P. sakaiensis* coding sequence reference using Salmon. Transcript-level counts were imported into R with tximport and analyzed with DESeq2 (Love et al. 2014). The condition factor was releveled so that maltose served as the reference level, making PET versus maltose the primary carbon-matched contrast. Differential expression was tested with a fold-change threshold of at least two-fold (lfcThreshold = 1) at a 5% false discovery rate (alpha = 0.05). Log2 fold changes were shrunk using apeglm (Zhu et al. 2019). Genes were annotated by parsing the coding sequence reference headers, and pathway genes were identified by searching protein descriptions for terms drawn from the source literature.

## Results

Growth on PET produced a large-scale transcriptional response. Principal component analysis separated the three PET replicates from the maltose and complex-medium samples along PC1, which captured 91% of the variance, while maltose and complex samples were largely indistinguishable (Figure 1). At the two-fold and 5% false discovery rate thresholds, 720 genes were differentially expressed between PET and maltose (389 up, 331 down).

The two enzymes responsible for the initial breakdown of PET were both significantly upregulated on PET (Figure 2). The gene encoding PETase showed the strongest induction of any pathway gene, increasing roughly 29-fold (log2 fold change 5.86, adjusted p = 1.5e-134). The gene encoding MHETase, which acts on the product of PETase, was also significantly upregulated, though more modestly, at roughly two-fold (log2 fold change 0.99, adjusted p = 7.1e-07). This pattern is consistent with induction of the two characterized PET-degrading enzymes of this organism (Yoshida et al. 2016) when plastic is present as a carbon source.

Genes in the downstream ethylene glycol and aromatic degradation pathways showed smaller and more mixed changes. A PQQ-dependent ethylene glycol dehydrogenase (PedE) and one subunit of protocatechuate 3,4-dioxygenase were modestly upregulated, whereas glycolate oxidase subunits and glyoxylate carboligase were downregulated (Figure 2). These smaller effects suggest that the strong, specific induction of PETase and MHETase is not matched by uniform upregulation of the entire downstream catabolic network under these conditions.

## Interpretation

The transcriptional response to PET is dominated by the front end of the degradation route. PETase, the secreted enzyme that attacks the plastic polymer directly, is induced far more strongly than any other pathway gene, and MHETase, which processes the intermediate PETase releases, is also significantly induced. Together this is consistent with the cell running the initial two-step breakdown of PET when plastic is the available carbon source (Yoshida et al. 2016).

The downstream pathways that process the resulting fragments, ethylene glycol and the aromatic terephthalate-derived compounds, show only small and mixed changes rather than coordinated induction. This suggests the response is selective rather than a wholesale activation of every gene in the pathway, though the smaller effect sizes and the limitations below mean this should be interpreted cautiously.

## Limitations

The PET cultures were sampled at approximately 27 hours versus approximately 9 hours for maltose (Novak and Gardner 2026), so the PET versus maltose contrast reflects a combination of carbon source and growth stage rather than carbon source alone. In addition, all conditions used 702 complex medium as a base, with maltose and PET added rather than substituted, so PET was not the sole carbon source. These design features of the source dataset limit interpretation of the broader transcriptional changes. The strong and specific induction of PETase, however, is unlikely to be explained by growth stage alone.

## Figures

**Figure 1.** Principal component analysis of all nine samples after variance-stabilizing transformation. `results/figures/pca_plot.png`

**Figure 2.** Log2 fold changes of significant PET-degradation pathway genes (PET versus maltose, adjusted p < 0.05). Values are log2 fold changes: each unit is a doubling, positive means higher expression on PET, negative means lower. `results/figures/pathway_genes.png`

## References

Love, M.I., Huber, W., and Anders, S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology. 15: 550. https://doi.org/10.1186/s13059-014-0550-8

Novak, J.K., and Gardner, J.G. 2026. Transcriptome data sets of *Piscinibacter sakaiensis* grown in 702 complex medium, maltose, and PET plastic. Microbiology Resource Announcements. 15: e01466-25. https://doi.org/10.1128/mra.01466-25

Yoshida, S., Hiraga, K., Takehana, T., Taniguchi, I., Yamaji, H., Maeda, Y., Toyohara, K., Miyamoto, K., Kimura, Y., and Oda, K. 2016. A bacterium that degrades and assimilates poly(ethylene terephthalate). Science. 351: 1196-1199. https://doi.org/10.1126/science.aad6359

Zhu, A., Ibrahim, J.G., and Love, M.I. 2019. Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. Bioinformatics. 35: 2084-2092. https://doi.org/10.1093/bioinformatics/bty895