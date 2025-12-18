import pandas as pd
import numpy as np
from pathlib import Path

HANDLE = "coroiusamy"
DATA_DIR = Path(f"data/work/{HANDLE}/lab10")
DATA_DIR.mkdir(parents=True, exist_ok=True)

np.random.seed(42) # Pentru rezultate reproductibile
n_samples = 50
n_snps = 20
n_genes = 30
samples = [f"Sample_{i:02d}" for i in range(1, n_samples + 1)]

# 1. Generare SNP Matrix (0, 1, 2)
snp_names = [f"rs{1000 + i}" for i in range(n_snps)]
snp_data = np.random.randint(0, 3, size=(n_samples, n_snps))
df_snp = pd.DataFrame(snp_data, index=samples, columns=snp_names)

# 2. Generare Expression Matrix (initializata cu zgomot random)
gene_names = [f"GENE_{i:02d}" for i in range(n_genes)]
exp_data = np.random.normal(0, 1, size=(n_samples, n_genes))
df_exp = pd.DataFrame(exp_data, index=samples, columns=gene_names)

# 3. Injectare Semnal (Corelatii artificiale)
# primele 5 gene sa fie puternic influentate de primele 5 SNP-uri
for i in range(5):
    snp_col = snp_names[i]
    gene_col = gene_names[i]
    
    # Formula: Expresie = 0.8 * SNP + 0.2 * Zgomot
    # Asta garanteaza o corelatie r > 0.7
    df_exp[gene_col] = (0.8 * df_snp[snp_col]) + np.random.normal(0, 0.2, n_samples)

# 4. Salvare fisiere
df_snp.to_csv(DATA_DIR / f"snp_matrix_{HANDLE}.csv")
df_exp.to_csv(DATA_DIR / f"expression_matrix_{HANDLE}.csv")
