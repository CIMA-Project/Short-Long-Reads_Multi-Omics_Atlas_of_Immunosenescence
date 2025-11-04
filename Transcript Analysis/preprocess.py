"""
This script preprocesses single-cell transcriptomic data using Scanpy and scCyclone.
"""
import anndata
import scanpy as sc
import scCyclone as scc
import numpy as np
import pandas as pd
import re
import scipy.sparse as sp

barcodes_path = './OUT.discovered_transcript_grouped_counts.barcodes.tsv'
features_path = './OUT.discovered_transcript_grouped_counts.features.tsv'
matrix_path = './OUT.discovered_transcript_grouped_counts.matrix.mtx'
sqanti_path = './OUT.novel.transcript_models_classification.txt'
gtf_file = "./OUT.novel.transcript_models_corrected.gtf"

def read_scCyclone_mtx(barcode_path,feature_path,matrix_path):
    sub_adata = anndata.read_mtx(matrix_path).T
    features = pd.read_csv(feature_path, header=None, sep='\t')
    sub_adata.var_names = anndata.utils.make_index_unique(pd.Index(features[0].values))
    sub_adata.obs_names = pd.read_csv(barcode_path, header=None)[0].values
    return sub_adata

adata = read_scCyclone_mtx(barcodes_path,features_path,matrix_path)
adata.var['isoform'] = adata.var.index
scc.tl.add_sqanti3(adata,sqanti_path)
if "predicted_NMD" in adata.var.columns:
    adata.var["predicted_NMD"] = adata.var["predicted_NMD"].astype(str)
adata = adata[:, adata.var['associated_gene'].dropna().index]
adata.var_names_make_unique()
sc.pp.filter_genes(adata,min_cells=3)
sc.pp.filter_cells(adata,min_genes=5)

transcript_id_to_name = {}
with open(gtf_file, 'r') as f:
    for line in f:
        match = re.search(r'transcript_id "(.+?)";.+gene_id "(.+?)";.+gene_name "(.+?)";', line)
        if match:
            transcript_id, _, gene_name = match.groups()
            transcript_id_to_name[transcript_id] = gene_name
adata.var['gene_name'] = pd.Series(
    adata.var_names, index=adata.var_names
).map(transcript_id_to_name)

mito_mask_index = ~adata.var_names.str.contains("chrM", case=False, na=False)
ribo_mask = ~adata.var["gene_name"].astype(str).fillna("").str.strip().str.upper().str.startswith(("RPS", "RPL"))
final_filter = mito_mask_index & ribo_mask
adata = adata[:, final_filter].copy()
keep_categories = ["full-splice_match", "novel_not_in_catalog", "incomplete-splice_match", "novel_in_catalog"]
adata = adata[:, adata.var["structural_category"].isin(keep_categories)]
adata.var['mt'] = adata.var['gene_name'].str.startswith('MT-')
adata.var['hb'] = adata.var['gene_name'].str.contains('^HB[^(P)]')
adata.var['rp'] = adata.var['gene_name'].str.match(r'^RP[SL][0-9]')
adata.var['ncRNA'] = adata.var['gene_name'].str.match(r'^[A-Z][A-Z][0-9].*\.[0-9]')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','rp','ncRNA','hb'], percent_top=None, log1p = False, inplace=True)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts','pct_counts_mt','pct_counts_rp','pct_counts_ncRNA','pct_counts_hb'], jitter=0.4, multi_panel=True)

samples = adata.obs['sample'].unique()
n_vars = adata.n_vars
gene_sample_count = np.zeros(n_vars, dtype=int)

for sample in samples:
    mask = (adata.obs['sample'] == sample).values
    X_sub = adata.X[mask, :]
    if sp.issparse(X_sub):
        expr = np.array((X_sub > 0).sum(axis=0)).ravel()
    else:
        expr = (X_sub > 0).sum(axis=0)
    gene_sample_count += (expr > 0).astype(int)
keep_idx = np.where(gene_sample_count >= 3)[0]
adata = adata[:, keep_idx].copy()

sr_adata = sc.read('./data_sr.h5ad')
mapping_df = sr_adata.obs[["celltype", "L2_celltype", "L3_celltype"]].copy()
common_cells = adata.obs_names.intersection(mapping_df.index)
adata = adata[common_cells].copy()
adata.obs[["celltype", "L2_celltype", "L3_celltype"]] = mapping_df.loc[common_cells]

adata.write('./data_lr.h5ad',compression='gzip')