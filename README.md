# Data Repo for Bactoridales paper at DFI

Title: Comprehensive analyses of a large human gut Bacteroidales culture collection reveal species and strain level diversity and evolution

## Folder structure

1. **code**: contains all R scripts for generating umaps, volcano plots and/or taxonomy barplots, etc. Scripts are organized by their order in the paper. Resultant plots are labeled with same number as the script. 
2. **results**: contains all unmodified graphs in PDF format. Graphs published along with the paper have been slightly rearranged and beautified in Adobe Illustrator.
3. **data**: contains raw, derived data or template files.
- `custom_bac71_phylogenomic-tree.txt`: Newick tree of representative isolates and their type strains based on curated single-copy core genes.
- `genome.size.csv`: list of representative isolates and their type strains.
- `illumina_typeStrain.prokka.rds`: parsed customized prokka output in R rds format. It only contains annotation information.
- `416_Illumina_only_nucmer_results.tsv`: parsed nucmer output of pairwise comparison of isolates of the same species.
- `new_output_from_CAZyme_hmmsearch.txt`: output from hmmsearch.
- `sppal.extended.rds`: species color palette in R rds format.
- `metab.quant.matrix.csv`: Production (positive values) or consumption (negative values) of SCFAs in mM of 111 isolates measured by quantitative metabolomics.
- `GMM_pathway_results.csv`: Completeness of metabolic pathways of each isolate, parsed result from [GOmixer](https://www.raeslab.org/omixer/).
- `metab.qual.matrix.log2.csv`: log2(fold change) of a panel of metabolites of 111 isolates measured by semi-quantitative metabolomics.
- `Final_ODs.xlsx`: Final OD600 readings of two biological replicates of each isolate measured for metabolomics.
- `donor.16S.csv`: 16S taxonomy profile for donors.
- `Anvio_geneCluster.matrix.csv`: Anvio gene cluster result in CSV format. 
- `pulpy.rds`: parsed result from [PULpy](https://github.com/WatsonLab/PULpy).
