# Config setup

## Configfile

Here 
- samples: Path to samples file 
	- format .tsv
	- required columns: sample
- counts: Path to count file
	- format: .tsv
	- format: first column gene names
	- column names == sample column from samples
-  gene_name_type:
	- system used for gene naming 
	- supports
		- ENSEMBL
		- HGNC
		- ENTREZ_ID
	- organisms:
		- Mus musculus
		- Homo sapiens
- dirs: Directories 
	- BASE_DATA_DIR: Directory where data is saved
	- BASE_ANALYSIS_DIR: Directory where results are saved

- pca:
	- labels: columns of sample sheet used for PCA plots (Used in snakemake report)
- diffexp: conditions for PCA
	- pval_threshold: padj threshold for DESeq2 test
	- LFC_threshold: absolute log2 fold change threshold 
	- contrasts: object containing the contrasts
		- <contrast_name>
			- <Comparison_one>
				- group_1
				- group_2
		- model: design formula definition for DESeq2 all variables need to be present in the sample .tsv file

- group_colors: colors to use in plots for groups
	- <contrast_name>
		- group: color
	- Color is accepted either as a name or hex code
		- "red" or "#ff0000"

- run_mitch: bool whether to run mitch on GSEA results
- run_carnival:
	- vanilla: bool Run vanilla carnival on DESeq2 results with perturbation targets (Custom gene support not added yet)
	- inverse: bool Run carnival on DESeq results without defined perturbation targets
	- sample: bool Run carnival on each sample separately
	- joint_rep: bool Generate joint HTML report for one contrast across all comparisons

- cplex_solver: Path to the cplex solver executable. 

## Sample sheet

THe sample sheet is required so that the workflow knows which condition is associated to each sample. 

Each sample should have an id which is identical in the count matrix and the first column of the sample sheet. After all contrast_names used in the config file need to have a column in the sample sheet. 

Have a look at data/test_data for example data used in the DESeq2 vignette. 