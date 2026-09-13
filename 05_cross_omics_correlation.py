import pandas as pd
import numpy as np
import pickle
from scipy import stats

TREATMENT_ORDER = ["C","H","D","S","HD","SD","SH","SHD"]
FDR_THRESHOLD = 0.05
LOG2FC_THRESHOLD = 1.0
R_THRESHOLD = 0.9
CORR_FDR_THRESHOLD = 0.05

def bh_fdr(pvals):
    p = np.asarray(pvals); n = len(p)
    order = np.argsort(p); ranked = p[order]
    fdr = ranked * n / (np.arange(n) + 1)
    fdr = np.minimum.accumulate(fdr[::-1])[::-1]
    fdr = np.clip(fdr, 0, 1)
    out = np.empty(n); out[order] = fdr
    return out

def main():
    with open("de_results.pkl", "rb") as f:
        de_results = pickle.load(f)
    with open("raw_means_filtered.pkl", "rb") as f:
        raw_means_f = pickle.load(f)
    with open("metab_means.pkl", "rb") as f:
        metab_means = pickle.load(f)

    shd = de_results["SHD"]
    shd_sig = shd[(shd.FDR < FDR_THRESHOLD) & (shd.log2FC.abs() >= LOG2FC_THRESHOLD)].index
    gene_log_means = np.log2(raw_means_f[TREATMENT_ORDER] + 1)
    gene_profiles = gene_log_means.loc[shd_sig]

    rows = []
    for metab in metab_means.columns:
        mvec = metab_means[metab].values
        if np.isnan(mvec).any():
            continue
        for gene in gene_profiles.index:
            gvec = gene_profiles.loc[gene].values
            r, p = stats.pearsonr(gvec, mvec)
            rows.append((gene, metab, r, p))

    corr_df = pd.DataFrame(rows, columns=["Gene","Metabolite","r","p"])
    corr_df["FDR"] = bh_fdr(corr_df["p"].values)
    strong = corr_df[(corr_df.r.abs() >= R_THRESHOLD) & (corr_df.FDR < CORR_FDR_THRESHOLD)]
    print(f"Significant pairs: {len(strong)}, unique genes: {strong.Gene.nunique()}")

    with open("corr_results.pkl", "wb") as f:
        pickle.dump({"corr_df": corr_df, "strong": strong}, f)

if __name__ == "__main__":
    main()
