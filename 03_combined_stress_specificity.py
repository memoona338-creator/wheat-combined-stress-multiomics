import pickle

FDR_THRESHOLD = 0.05
LOG2FC_THRESHOLD = 1.0

def significant_genes(de_results, treatment, fdr=FDR_THRESHOLD, fc=LOG2FC_THRESHOLD):
    r = de_results[treatment]
    return set(r[(r.FDR < fdr) & (r.log2FC.abs() >= fc)].index)

def main():
    with open("de_results.pkl", "rb") as f:
        de_results = pickle.load(f)

    sig = {t: significant_genes(de_results, t) for t in de_results}
    for t, genes in sig.items():
        print(f"{t}: {len(genes)} significant genes")

    shd_unique = sig["SHD"] - sig["H"] - sig["D"] - sig["S"] - sig["HD"] - sig["SD"] - sig["SH"]
    pct = 100 * len(shd_unique) / len(sig["SHD"]) if sig["SHD"] else 0
    print(f"SHD-unique genes: {len(shd_unique)} ({pct:.1f}% of {len(sig['SHD'])})")

    with open("sig_sets.pkl", "wb") as f:
        pickle.dump({"sig": sig, "shd_unique": shd_unique}, f)

if __name__ == "__main__":
    main()
