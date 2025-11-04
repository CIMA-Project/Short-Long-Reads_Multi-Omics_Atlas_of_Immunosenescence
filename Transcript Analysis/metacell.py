"""
This script processes single-cell RNA-seq data to assign metacells using SEACells and summarizes gene expression by metacell.
"""
import numpy as np
import pandas as pd
import anndata as ad
from tqdm import tqdm
from sklearn.preprocessing import normalize
from sklearn.decomposition import TruncatedSVD         
import scipy.sparse
from SEACells import core
import scanpy as sc
import psutil, os, sys, time, gc

PROC = psutil.Process(os.getpid())

def mem(tag):
    """Print current RSS(GB)"""
    rss = PROC.memory_info().rss / 1024**3
    print(f"[{time.strftime('%H:%M:%S')}] {tag} | RSS = {rss:,.1f} GB", file=sys.stderr)
    sys.stderr.flush()
    return rss

os.makedirs("metacell", exist_ok=True)

def tfidf_safe(X):
    """TF-IDF with guards against 0/0→NaN."""
    idf = X.shape[0] / np.maximum(1, X.sum(axis=0))
    if scipy.sparse.issparse(X):
        tf = X.multiply(1 / np.maximum(1, X.sum(axis=1)))
        X = tf.multiply(idf)
        X.data = np.nan_to_num(X.data, copy=False)
    else:
        tf = X / np.maximum(1, X.sum(axis=1, keepdims=True))
        X = tf * idf
        np.nan_to_num(X, copy=False)
    return X

def assign_metacells(
        ad, samples="sample", celltypes="L3_celltype",
        svd_dim=50, cells_per_meta=5):

    for sample in tqdm(ad.obs[samples].unique(), desc=samples):
        tmp_sample = ad[ad.obs[samples] == sample]
        for ctype in tqdm(tmp_sample.obs[celltypes].unique(),
                          desc=celltypes, leave=False):
            sub = tmp_sample[tmp_sample.obs[celltypes] == ctype].copy()
            print(f"\n### sample={sample} ctype={ctype} n={sub.n_obs}", file=sys.stderr)

            if sub.n_obs < 20:
                sub.obs["SEACell_ID"] = f"{sample}-{ctype}-SEACell-0"
                ad.obs.loc[sub.obs_names, "SEACell_ID"] = sub.obs["SEACell_ID"]
                continue

            X = tfidf_safe(sub.X)
            X = normalize(X, norm="l1", axis=1, copy=False)
            # Multiply by 1e4 to scale normalized values before log transformation, following common scRNA-seq preprocessing practice
            X = np.log1p(X * 1e4)

            svd = TruncatedSVD(n_components=svd_dim, random_state=666)
            X_lsi = svd.fit_transform(X)
            X_lsi -= X_lsi.mean(1, keepdims=True)
            X_lsi /= np.maximum(1e-6, X_lsi.std(1, ddof=1, keepdims=True))
            sub.obsm["X_lsi"] = X_lsi

            n_SEACells = int(np.floor(sub.n_obs / cells_per_meta)) + 1

            nw_eigs = max(1, min(4, X_lsi.shape[1]))

            if sub.n_obs <= n_SEACells or sub.n_obs <= nw_eigs:
                sub.obs["SEACell_ID"] = f"{sample}-{ctype}-SEACell-0"
                ad.obs.loc[sub.obs_names, "SEACell_ID"] = sub.obs["SEACell_ID"]
                continue

            model = core.SEACells(
                sub,
                build_kernel_on="X_lsi",
                n_SEACells=n_SEACells,
                n_waypoint_eigs=nw_eigs,
                convergence_epsilon=1e-5,
                use_gpu=False,
            )
            model.construct_kernel_matrix()
            model.initialize_archetypes()
            try:
                model.fit(min_iter=10, max_iter=500)

            except RuntimeWarning:
                model = core.SEACells(
                    sub,
                    build_kernel_on="X_lsi",
                    n_SEACells=max(1, int(n_SEACells * 0.5)),
                    n_waypoint_eigs=nw_eigs,
                    convergence_epsilon=5e-5,
                    use_gpu=False,
                )
                model.construct_kernel_matrix()
                try:
                    model.initialize_archetypes()
                    model.fit(min_iter=10, max_iter=1000)
                except (RuntimeWarning, IndexError):
                    sub.obs["SEACell_ID"] = f"{sample}-{ctype}-SEACell-0"
                    ad.obs.loc[sub.obs_names, "SEACell_ID"] = sub.obs["SEACell_ID"]
                    del sub, X, X_lsi, svd, model
                    continue

            sub.obs["SEACell_ID"] = (
                f"{sample}-{ctype}-" + sub.obs["SEACell"].astype(str)
            )
            ad.obs.loc[sub.obs_names, "SEACell_ID"] = sub.obs["SEACell_ID"]

            del sub, X, X_lsi, svd, model
            gc.collect()

    return ad

def summarize_by_SEACell(ad_in, label="SEACell"):
    from scipy.sparse import csr_matrix
    metas = ad_in.obs[label].unique()
    df = pd.DataFrame(0.0, index=metas, columns=ad_in.var_names)
    for mc in tqdm(df.index, desc="Summarizing"):
        cells = ad_in.obs_names[ad_in.obs[label] == mc]
        df.loc[mc, :] = np.ravel(ad_in[cells, :].X.sum(axis=0))
    meta_adata = ad.AnnData(csr_matrix(df.values))
    meta_adata.obs_names = df.index.astype(str)
    meta_adata.var_names = df.columns
    return meta_adata

lr = sc.read_h5ad("./data_lr.h5ad")
lr = assign_metacells(lr)
lr.obs["SEACell"] = lr.obs["SEACell_ID"].astype("category")
lr.obs[["SEACell"]].to_csv("./metacell_assignment.csv")
lr_meta = summarize_by_SEACell(lr)
lr_meta.write("./metacell.h5ad")
