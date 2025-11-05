## ------------------------------
##tcr preprocessing
## ------------------------------
adata_t = sc.read_h5ad('ALSM/RNA/ALL_T_L3_annotation.h5ad')
path = 'ALSM/TCR/data_with_flags/'
samples = [s for s in os.listdir(path)][1:]

def read_all_tcr(sample):
    sample_path = os.path.join(path, sample)
    adata_merge = ir.io.read_10x_vdj(sample_path)
    adata_merge['reads'] = adata_merge['umis']
    sample_id = sample.replace('_TCR_filtered_contig_annotations.csv', '')
    adata_merge.obs_names = sample_id + '_' + adata_merge.obs_names
    print(f"{sample} -- success! {adata_merge.shape}")
    return adata_merge

adata_list = []
for sample in samples:
    adata_list.append(read_all_tcr(sample))
# 合并所有 AnnData
adata_tcr = anndata.concat(adata_list, join="outer")
print("合并完成：", adata_tcr.shape)
##muon 容器
mdata_tcr = mu.MuData({"gex": adata_t, "airr": adata_tcr})
ir.pp.index_chains(mdata_tcr)
ir.tl.chain_qc(mdata_tcr)

print(mdata_tcr['airr'].obs['chain_pairing'].value_counts(),
(mdata_tcr["airr"].obs["chain_pairing"]
    .value_counts(normalize=True)          
    .mul(100)                             
    .round(6)))
##质控
mu.pp.filter_obs(mdata_tcr['airr'], "chain_pairing", lambda x: ~np.isin(x, ["orphan VDJ","orphan VJ","ambiguous"]))
##TCR clonotype_network L3_celltype
ir.pp.ir_dist(mdata_tcr,sequence="aa")
ir.tl.define_clonotypes(mdata_tcr, receptor_arms="all",dual_ir="primary_only")
ir.tl.clonotype_network(mdata_tcr, min_cells=3)
L3_color_map_t = {
    'CD4 Naïve': "#015696", 'CD4 Regulatory': "#89CAEA", 'CD4 CTL': "#BBE6FA",
    'CD4 Tcm': "#0B75B3", 'CD4 Tem': "#4596CD",   
    'NKT': "#87D649", 'MAIT': "#47A183", 'γδT': "#24858E", 'Cycling T': "#77C3BC",
    'CD8 CTL': "#C3E023", 'CD8 Naïve': "#8DCA9D", 'CD8 Tem': "#2BB17E", 
    'CD8 Tcm': "#50C56A"
}
clonotype_network_clonotype_L3 = ir.pl.clonotype_network(mdata_tcr, 
                            color="gex:L3_celltype",
                            palette=L3_color_map_t,
                            base_size=20, 
                            title="TCR Clonotype Network by L3_celltype(more than 3 cells)",
                            label_fontsize=9, 
                            panel_size=(12, 8))
clonotype_network_clonotype_L3.figure.savefig('TCR/pdf/clonotype_network_L3.pdf', 
                                     bbox_inches='tight', 
                                     dpi=300)
##TCR clonotype_network Age                                    
group_colors = {
    "OLD": "#d88f91",
    "YOUNG": "#80bcc8"}
clonotype_network_clonotype = ir.pl.clonotype_network(mdata, 
                            color="gex:group",
                            palette=group_colors,
                            base_size=20, 
                            title="TCR Clonotype Network by Age(more than 3 cells)",
                            label_fontsize=9, 
                            panel_size=(12, 8))
clonotype_network_clonotype.figure.savefig('TCR/pdf/clonotype_network_Age.pdf', 
                                     bbox_inches='tight', 
                                     dpi=300)

## ------------------------------
##bcr preprcessing
## ------------------------------
adata_b =sc.read_h5ad('ASLM/RNA/B_L3_annotation_cellID.h5ad')
path = 'ALSM/BCR/data/'
samples = [s for s in os.listdir(path)]

def read_all_bcr(sample):
    sample_path = os.path.join(path, sample)
    adata_merge = ir.io.read_10x_vdj(sample_path)
    sample_id = sample.replace('_BCR_filtered_contig_annotations.csv', '')
    adata_merge['reads'] = adata_merge['umis']
    adata_merge.obs_names = sample_id + '_' + adata_merge.obs_names
    print(f"{sample} -- success! {adata_merge.shape}")
    return adata_merge
adata_list = []
for sample in samples:
    adata_list.append(read_all_bcr(sample))

# 合并所有 AnnData
adata_bcr = anndata.concat(adata_list, join="outer")
print("合并完成：", adata_bcr.shape)
mdata_bcr = mu.MuData({"gex": adata_b, "airr": adata_bcr})
ir.pp.index_chains(mdata_bcr)
ir.tl.chain_qc(mdata_bcr)

print(mdata_bcr['airr'].obs['chain_pairing'].value_counts(),
(mdata_bcr["airr"].obs["chain_pairing"]
    .value_counts(normalize=True)          
    .mul(100)                             
    .round(6)))
mu.pp.filter_obs(mdata_bcr['airr'], "chain_pairing", lambda x: ~np.isin(x, ["orphan VDJ","orphan VJ"]))
## convert_to_dandelion
vdjx_1 = ddl.from_scirpy(adata_bcr)
vdjx_1.data["v_call"].replace("", np.nan, inplace=True)
vdjx_1.data.dropna(subset=["v_call"], inplace=True)

vdjx_1.data["j_call"].replace("", np.nan, inplace=True)
vdjx_1.data.dropna(subset=["j_call"], inplace=True)

vdjx_1.data["junction_aa"].replace("", np.nan, inplace=True)
vdjx_1.data.dropna(subset=["junction_aa"], inplace=True)

vdjx_1.data["junction_length"] = [len(a) for a in vdjx_1.data["junction_aa"]]
ddl.pp.calculate_threshold(vdjx_1, model="hh_s5f")

ir.pp.ir_dist(
    mdata_bcr,
    metric="normalized_hamming",
    sequence="nt",
    cutoff=vdjx.threshold*100,)

ir.tl.define_clonotype_clusters(
    mdata_bcr,
    sequence="nt",
    metric="normalized_hamming",
    receptor_arms="all",
    dual_ir="primary_only",
    same_v_gene=True,
    same_j_gene=True,
    partitions="fastgreedy",
    key_added="clone_id_similarity",
)

ir.tl.clonotype_network(
    mdata_airr, sequence="nt", metric="normalized_hamming", min_cells=3, clonotype_key="clone_id_similarity"
)
## cloneNetwork age
group_colors = {
    "OLD": "#d88f91",
    "YOUNG": "#80bcc8"}
BCR_clonotype_network_clonotype_Age = ir.pl.clonotype_network(
    mdata_airr,
    color="gex:group",
    palette=group_colors,
    title="Clonotype Network BCR by Age(more than 3 cells)",
    label_fontsize=9,
    panel_size=(12, 8),
    base_size=20,
)
BCR_clonotype_network_clonotype_Age.figure.savefig('BCR/pdf/clonotype_network_Age_group.pdf', 
                                     bbox_inches='tight', 
                                     dpi=300)
## cloneNetwork L3
L3_color_map_b = {
    'Swit. Bm': "#D48CB1", 'Unswit. Bm': "#8085C0", 'Naïve B': "#D4A7C4", 
    'Plasma': "#DBC6D9", 'Plsamablast': "#E6E6F0",
}
BCR_clonotype_network_clonotype_L3 = ir.pl.clonotype_network(
    mdata_airr,
    color="gex:L3_celltype",
    palette=L3_color_map_b,
    title="Clonotype Network BCR by L3_celltype(more than 3 cells)",
    label_fontsize=9,
    panel_size=(12, 8),
    base_size=20,
)
BCR_clonotype_network_clonotype_L3.figure.savefig('BCR/pdf/clonotype_network_L3.pdf', 
                                     bbox_inches='tight', 
                                     dpi=300)



##TRA CDR3 Length(aa) by L3_celltype
cdr3_length_L3 = ir.pl.spectratype(
    mdata,
    color="gex:L3_celltype",
    viztype="curve",
    chain='VJ_1',
    palette=L3_color_map_t,
    curve_layout="shifted",
    fig_kws={"dpi": 300},
    kde_kws={"kde_norm": False},
    title="cdr3 Length by L3_celltype"
)
cdr3_length_L3.figure.savefig('TCR/pdf/TRA_CDR3aaLength_L3.pdf', 
                                     bbox_inches='tight', 
                                     dpi=300)


