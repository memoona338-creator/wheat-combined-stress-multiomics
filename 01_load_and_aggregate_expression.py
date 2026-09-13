import pandas as pd
import pickle
import re

INPUT_FILE = "Data_S2_TPM.xlsx"
SHEET_NAME = "Expression Data_TPM"
HEADER_ROW = 6

def gene_id_from_transcript(transcript_id):
    t = str(transcript_id)
    t = re.sub(r"LC\.\d+$", "", t)
    t = re.sub(r"\.\d+$", "", t)
    return t

def main():
    tpm = pd.read_excel(INPUT_FILE, sheet_name=SHEET_NAME, header=HEADER_ROW)
    tpm.rename(columns={tpm.columns[0]: "Transcript_ID"}, inplace=True)
    print(f"Loaded transcript-level TPM matrix: {tpm.shape}")

    tpm["Gene_ID"] = tpm["Transcript_ID"].apply(gene_id_from_transcript)
    sample_cols = [c for c in tpm.columns if c not in ("Transcript_ID", "Gene_ID")]
    gene_tpm = tpm.groupby("Gene_ID")[sample_cols].sum()
    print(f"Aggregated to gene level: {gene_tpm.shape}")

    with open("gene_tpm.pkl", "wb") as f:
        pickle.dump(gene_tpm, f)

if __name__ == "__main__":
    main()
