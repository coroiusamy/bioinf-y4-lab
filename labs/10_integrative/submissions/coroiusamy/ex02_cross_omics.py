"""
Exercise 10.2 — Identify top SNP–Gene correlations

TODO:
- încărcați matricea integrată multi-omics
- împărțiți rândurile în SNPs vs gene (după indice sau după nume)
- calculați corelații între fiecare SNP și fiecare genă
- filtrați |r| > 0.5
- exportați snp_gene_pairs_<handle>.csv
"""

from pathlib import Path
import pandas as pd

HANDLE = "coroiusamy"
JOINT_CSV = Path(f"labs/10_integrative/submissions/{HANDLE}/multiomics_concat_{HANDLE}.csv")

OUT_CSV = Path(f"labs/10_integrative/submissions/{HANDLE}/snp_gene_pairs_{HANDLE}.csv")

# TODO: load joint matrix, compute correlations, export

# 1. Incarcare matrice integrata
df_joint = pd.read_csv(JOINT_CSV, index_col=0)

# 2. Impartire coloane in SNP si Gene
snp_cols = [c for c in df_joint.columns if c.startswith('rs')]
gene_cols = [c for c in df_joint.columns if not c.startswith('rs')]

# 3. Calcul corelatii cross-omics
correlations = []
for snp in snp_cols:
    for gene in gene_cols:
        r = df_joint[snp].corr(df_joint[gene])
        
        # 4. Filtrare prag |r| > 0.5
        if abs(r) > 0.5:
            correlations.append({
                'SNP': snp,
                'Gene': gene,
                'Correlation': r
            })

# 4. Export rezultate
df_res = pd.DataFrame(correlations)
if not df_res.empty:
    df_res = df_res.sort_values(by='Correlation', key=abs, ascending=False)
df_res.to_csv(OUT_CSV, index=False)