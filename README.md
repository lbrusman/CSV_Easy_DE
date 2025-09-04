# CSV_Easy_DE
This is Liza's pipeline for running DESeq2 and RUVseq from a csv input

# Here is how you run this pipeline:
## Step 1:

Clone this repo to the directory where you want to do your analysis.

## Step 2:
Make sure your files/directories are in the following format:

- `metadata`: Must be a .tsv or .txt file where each row represents one sample.
- `indir`: Must be a directory with all pseudobulk matrix files (with **only** files you want to run).
- `pseudobulk matrices`: Must be .tsv or .txt files where each column represents one sample and each row represents one feature (gene). Must all have filenames in the format `<celltype>_filepattern.tsv`. File pattern must be the same across all files.

## Step 3:
**Important note:** metadata column names must be free of underscores `("_")`. Please change to periods `(".")` or eliminate if necessary before running this pipeline.

Modify [`input_csv.csv`](input_csv/input_csv.csv) to be compatible with your metadata file.
Here are descriptions of each column:

- `contrast_id`: what you would like to name the contrast. Generally used as output filename prefix.
- `contrast_var`: the name of the column you would like to perform the contrast on.
- `donor_col`: the name of the column with sample (donor) names.
- `control_grp`: if this is a binary contrast, whichever group is your control group (i.e. if expression is higher in control_grp, log2FoldChange will be negative).
- `experimental_grp`: if this is a binary contrast, whichever group is your experimental group (i.e. if expression is higher in experimetal group, log2FoldChange will be positive).
- `subset_on1`: would you like to subset data on a variable before running differential expression? Put that column name here.
- `subset_on2`: NOT YET IMPLEMENTED
- `subset1`: the groups you would like to *keep* from the column in `subset_on1`. The groups must be underscore-separated.
- `subset2`: NOT YET IMPLEMENTED
- `covariates`: names of columns (covariates) you would like to include in the differential expression formula. If more than one, must be underscore-separated.
- `correlate_vars`: names of columns you would like to correlate with latent variables in a correlation matrix.
- `scale`: do you want to scale numeric variables? Must be TRUE or FALSE.
- `file_pattern`: the file ending of all of your pseudobulk matrix files. This should include **everything** other than the cell type.
- `assay`: must be either "RNA" or "ATAC" (no quotes). Pipeline will pick parameters accordingly.

## Step 4:

Modify [`submit_jobs.sh`](src/submit_jobs.sh) to point to files and directories you want. 

**Note:** `input_csv.csv` path must be changed at the beginning AND end of the file.

## Step 5:

Modify [`run_jobs.sh`](src/run_jobs.sh) to point to your conda environment (with DESeq2 and RUVseq installed) and change all the SLURM info to your info.

## Step 6:

To run the pipeline:

`cd src`

`bash submit_jobs.sh`

This will submit a separate job for each contrast (row in your input csv).

