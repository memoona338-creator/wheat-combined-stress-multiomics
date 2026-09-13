import pandas as pd
import numpy as np
import pickle
from scipy import stats

TREATMENT_GROUPS = {
    "C": ["C1","C2","C4","C8"], "H": ["H1","H2","H3","H4"],
    "D": ["D1","D4","D6","D8"], "S": ["S1","S3","S5","S8"],
    "HD": ["HD3","HD4","HD7","HD8"], "SD": ["SD1","SD5","SD6","SD8"],
    "SH": ["SH2","SH3","SH4","SH8"], "SHD": ["SHD4","SHD5","SHD6","SHD7"],
}
CONTROL = "C"
FDR_THRESHOLD = 0.05
LOG2FC_THRESHOLD = 1.0

def bh_fdr(pvals):
    p = np.asarray(pvals); n = len(p)
    order = np.argsort(p); ranked = p[order]
    fdr = ranked * n / (np.arange(n) + 1)
    fdr = np.minimum.accumulate(fdr[::-1])[::-1]
    fdr = np.clip(fdr, 0, 1)
    out = np.empty(n); out[order] = fdr
    return out

def differential_expression(log_tpm, means, treat, ctrl=CONTROL):
    a = log_tpm[TREATMENT_GROUPS[treat]].values
    b = log_tpm[TREATMENT_GROUPS[ctrl]].values
    t_stat, p = stats.ttest_ind(a, b, axis=1, equal_var=False)
    p = np.nan_to_num(p, nan=1.0)
    fdr = bh_fdr(p)
    log2_means = np.log2(means + 1)
    l2fc = log2_means[treat] - log2_means[ctrl]
    return pd.DataFrame({"log2FC": l2fc.values, "pval": p, "FDR": fdr}, index=log_tpm.index)

def main():
    with open("gene_tpm.pkl", "rb") as f:
        gene_tpm = pickle.load(f)

    log_tpm = np.log2(gene_tpm + 1)
    raw_means = pd.DataFrame({g: gene_tpm[cols].mean(axis=1) for g, cols in TREATMENT_GROUPS.items()})
    expressed = raw_means.max(axis=1) > 1.0
    log_tpm_f = log_tpm[expressed]
    raw_means_f = raw_means[expressed]
    print(f"Expressed genes retained: {expressed.sum()}")

    results = {}
    for treat in [t for t in TREATMENT_GROUPS if t != CONTROL]:
        results[treat] = differential_expression(log_tpm_f, raw_means_f, treat)
        sig = results[treat][(results[treat].FDR < FDR_THRESHOLD) & (results[treat].log2FC.abs() >= LOG2FC_THRESHOLD)]
        print(f"{treat}: {len(sig)} DEGs  up={(sig.log2FC>0).sum()} down={(sig.log2FC<0).sum()}")

    with open("de_results.pkl", "wb") as f:
        pickle.dump(results, f)
    with open("raw_means_filtered.pkl", "wb") as f:
        pickle.dump(raw_means_f, f)

if __name__ == "__main__":
    main()
