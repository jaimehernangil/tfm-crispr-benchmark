import pandas as pd
from pathlib import Path

input_dir = Path(".")
output_file = "CHOPCHOP_top3_dataset.csv"

all_results = []

for file in input_dir.glob("*.txt"):
    gene = file.stem.replace("_tabla", "").replace("_CHOPCHOP", "")

    df = pd.read_csv(file, sep="\t")

    top3 = df.head(3).copy()
    top3.insert(0, "Gene", gene)

    top3["sgRNA_20nt"] = top3["Target sequence"].str[:-3]
    top3["PAM"] = top3["Target sequence"].str[-3:]

    all_results.append(top3)

final_df = pd.concat(all_results, ignore_index=True)

cols = [
    "Gene",
    "Rank",
    "sgRNA_20nt",
    "PAM",
    "Target sequence",
    "Genomic location",
    "Strand",
    "GC content (%)",
    "Self-complementarity",
    "MM0",
    "MM1",
    "MM2",
    "MM3",
    "Efficiency"
]

final_df = final_df[cols]

final_df.to_csv(output_file, index=False)

print("Dataset generado correctamente")
