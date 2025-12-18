"""
Exercise 10 — PCA Single-Omics vs Joint

TODO:
- încărcați SNP și Expression
- normalizați fiecare strat (z-score)
- rulați PCA pe:
    1) strat SNP
    2) strat Expression
    3) strat Joint (concat)
- generați 3 figuri PNG
- comparați vizual distribuția probelor
"""

from pathlib import Path
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

HANDLE = "coroiusamy"

SNP_CSV = Path(f"data/work/{HANDLE}/lab10/snp_matrix_{HANDLE}.csv")
EXP_CSV = Path(f"data/work/{HANDLE}/lab10/expression_matrix_{HANDLE}.csv")

OUT_DIR = Path(f"labs/10_integrative/submissions/{HANDLE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# TODO: load, align, normalize, PCA, export figures

# 1. Incarcare date
df_snp = pd.read_csv(SNP_CSV, index_col=0)
df_exp = pd.read_csv(EXP_CSV, index_col=0)

# 2. Aliniere probe comune
common_samples = df_snp.index.intersection(df_exp.index)
df_snp = df_snp.loc[common_samples]
df_exp = df_exp.loc[common_samples]

# 3. Normalizare z-score
scaler = StandardScaler()
snp_scaled = scaler.fit_transform(df_snp)
exp_scaled = scaler.fit_transform(df_exp)

# 4. Creare matrice Joint si salvare sub forma de CSV
df_joint = pd.concat([
    pd.DataFrame(snp_scaled, index=common_samples, columns=df_snp.columns),
    pd.DataFrame(exp_scaled, index=common_samples, columns=df_exp.columns)
], axis=1)
df_joint.to_csv(OUT_DIR / f"multiomics_concat_{HANDLE}.csv")

# 5. Functie pentru PCA si salvare PNG
def run_and_save_pca(data, title, filename):
    pca = PCA(n_components=2)
    coords = pca.fit_transform(data)
    plt.figure(figsize=(8, 6))
    plt.scatter(coords[:, 0], coords[:, 1], alpha=0.7)
    plt.title(title)
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.2%})")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.2%})")
    plt.grid(True)
    plt.savefig(OUT_DIR / filename)
    plt.close()

# Executie PCA pentru cele 3 variante
run_and_save_pca(snp_scaled, "PCA SNP Layer", "pca_snp.png")
run_and_save_pca(exp_scaled, "PCA Expression Layer", "pca_expr.png")
run_and_save_pca(df_joint, "PCA Joint Layer", f"pca_joint_{HANDLE}.png")