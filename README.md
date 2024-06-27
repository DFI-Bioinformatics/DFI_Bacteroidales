# Data Repo for Bactoridales paper at DFI

Title: Comprehensive analyses of a large human gut Bacteroidales culture collection reveal species and strain level diversity and evolution

## Folder structure

1. **code**: contains all R scripts for generating umaps, volcano plots and/or taxonomy barplots, etc. Scripts are organized by their order in the paper. Resultant plots are labeled with same number as the script. 
2. **results**: contains all unmodified graphs in PDF format. Graphs published along with the paper have been slightly rearranged and beautified in Adobe Illustrator.
3. **data**: contains raw, derived data or template files.
- `genome.size.csv`: list of representative isolates and their type strain.
- `illumina_typeStrain.prokka.rds`: parsed customized prokka output in R rds format. It only contains annotation information.
- `new_output_from_CAZyme_hmmsearch.txt`: output from hmmsearch.
- `sppal.extended.rds`: species color palette in R rds format.
- `donor.16S.csv`: 16S taxonomy profile for donors.
- `Anvio_geneCluster.matrix.csv`: Anvio gene cluster result in CSV format. 
- `pulpy.rds`: parsed result from [PULpy](https://github.com/WatsonLab/PULpy).
