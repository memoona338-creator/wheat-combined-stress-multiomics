import pandas as pd
import pickle

INPUT_FILE = "Data_S1_phenotypic_metabolic.xlsx"
HEADER_ROW = 4
VALID_TREATMENT_CODES = {"C","H","D","S","HD","SD","SH","SHD"}
TREATMENT_ORDER = ["C","H","D","S","HD","SD","SH","SHD"]

def clean_raw_rows(df, treat_col="Treatment"):
    df = df.copy()
    df[treat_col] = df[treat_col].astype(str).str.strip()
    return df[df[treat_col].isin(VALID_TREATMENT_CODES)]

def main():
    xl = pd.ExcelFile(INPUT_FILE)

    aa = clean_raw_rows(xl.parse("Amino_acids_data", header=HEADER_ROW))
    suc = clean_raw_rows(xl.parse("Sucrose_data", header=HEADER_ROW))[["Treatment","Sucrose (ug/mg)"]]
    oa = clean_raw_rows(xl.parse("Organic_acids", header=HEADER_ROW))[["Treatment","Total_organic_acids (ug/mg)"]]

    aa_means = aa.groupby("Treatment").mean(numeric_only=True).reindex(TREATMENT_ORDER)
    suc_means = suc.groupby("Treatment").mean(numeric_only=True).reindex(TREATMENT_ORDER)
    oa_means = oa.groupby("Treatment").mean(numeric_only=True).reindex(TREATMENT_ORDER)

    metab_means = pd.concat([aa_means, suc_means, oa_means], axis=1)
    print(f"Metabolite table: {metab_means.shape}")

    with open("metab_means.pkl", "wb") as f:
        pickle.dump(metab_means, f)

if __name__ == "__main__":
    main()
