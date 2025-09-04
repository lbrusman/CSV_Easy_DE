#make sure file can be read in properly
input_vars="../input_csv/input_csv.csv"

dos2unix "$input_vars" #change to unix

sed -i -e '$a\' "$input_vars" #add blank line at end of csv if not already there

#now run script for each contrast
while IFS=, read -r contrast_id contrast_var donor_col control_grp experimental_grp subset_on1 subset_on2 subset1 subset2 covariates correlate_vars scale file_pattern assay
do

echo $contrast_id
echo $contrast_var

export contrast_id
export contrast_var
export donor_col
export control_grp
export experimental_grp
export subset_on1
export subset_on2
export subset1
export subset2
export covariates
export correlate_vars
export scale
export file_pattern
export assay

indir=/tscc/projects/ps-gaultonlab/welison/mega.panc/results/pseudobulk.matrices/Pancreas/RNA/
export indir

outdir=/tscc/projects/ps-gaultonlab/lebrusman/mega_panc/RUVseq/outputs_filt_W_padj_max_k10_250825_old_contrasts/
mkdir -p $outdir
export outdir

metadata=/tscc/projects/ps-gaultonlab/welison/mega.panc/results/pseudobulk.matrices/Pancreas/RNA/Pancreas_RNA_DESeq_Meta_v1.tsv
export metadata


sbatch run_jobs.sh
done < <(tail -n +2 ../input_csv/input_csv.csv)
