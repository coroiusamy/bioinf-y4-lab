"""
Exercise 9.2 — Disease Proximity and Drug Ranking

Scop:
- să calculați distanța medie dintre fiecare medicament și un set de gene asociate unei boli
- să ordonați medicamentele în funcție de proximitate (network-based prioritization)

TODO-uri principale:
- încărcați graful bipartit drug–gene (din exercițiul 9.1) sau reconstruiți-l
- încărcați lista de disease genes
- pentru fiecare medicament, calculați distanța minimă / medie până la genele bolii
- exportați un tabel cu medicamente și scorul lor de proximitate
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Set, List, Tuple
import networkx as nx
import pandas as pd
import numpy as np
import pickle

# --------------------------
# Config
# --------------------------
HANDLE = "coroiusamy"

# Input
GRAPH_DRUG_GENE = Path(f"labs/09_repurposing/submissions/{HANDLE}/network_drug_gene_{HANDLE}.pkl")
DRUG_GENE_CSV = Path(f"data/work/{HANDLE}/lab09/drug_gene_{HANDLE}.csv") 
DISEASE_GENES_TXT = Path(f"data/work/{HANDLE}/lab09/disease_genes_{HANDLE}.txt")

# Output
OUT_DIR = Path(f"labs/09_repurposing/submissions/{HANDLE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_DRUG_PRIORITY = OUT_DIR / f"drug_priority_{HANDLE}.csv"

# --------------------------
# Utils
# --------------------------
def ensure_exists(path: Path) -> None:
    """
    TODO:
    - verificați că fișierul există
    - dacă nu, ridicați FileNotFoundError
    """
    if not path.exists():
        raise FileNotFoundError(f"Fișierul lipsește: {path}")


def load_bipartite_graph_or_build() -> nx.Graph:
    """
    TODO:
    - dacă GRAPH_DRUG_GENE există, încărcați-l direct
    - altfel, reconstruiți graful plecând de la DRUG_GENE_CSV
    """
    if GRAPH_DRUG_GENE.exists():
        with open(GRAPH_DRUG_GENE, 'rb') as f:
            return pickle.load(f)
    
    df = pd.read_csv(DRUG_GENE_CSV).dropna().astype(str)
    B = nx.Graph()
    for _, row in df.iterrows():
        B.add_node(row['drug'], bipartite="drug")
        B.add_node(row['gene'], bipartite="gene")
        B.add_edge(row['drug'], row['gene'])
    return B


def load_disease_genes(path: Path) -> Set[str]:
    """
    TODO:
    - încărcați fișierul text cu gene (una pe linie)
    - returnați un set de gene (string)
    """
    with open(path, 'r') as f:
        return {line.strip() for line in f if line.strip()}


def get_drug_nodes(B: nx.Graph) -> List[str]:
    """
    TODO:
    - extrageți lista nodurilor de tip 'drug'
    - presupunem atributul bipartite="drug"
    """
    return [n for n, d in B.nodes(data=True) if d.get("bipartite") == "drug"]


def compute_drug_disease_distance(
    B: nx.Graph,
    drug: str,
    disease_genes: Set[str],
    mode: str = "mean",
    max_dist: int = 8,
) -> float:
    """
    TODO:
    - calculați distanța shortest_path_length până la genele bolii
    """
    all_nodes = list(B.nodes)
    target_nodes = []
    # Matching inteligent pentru numele lungi din Open Targets
    for dg in disease_genes:
        matches = [node for node in all_nodes if dg.upper() in node.upper()]
        target_nodes.extend(matches)
    
    target_nodes = list(set(target_nodes))
    if not target_nodes: return float(max_dist)
    
    dists = []
    for t in target_nodes:
        try:
            dists.append(nx.shortest_path_length(B, drug, t))
        except nx.NetworkXNoPath:
            dists.append(float(max_dist))
    
    return np.mean(dists) if mode == "mean" else np.min(dists)


def rank_drugs_by_proximity(
    B: nx.Graph,
    disease_genes: Set[str],
    mode: str = "mean",
) -> pd.DataFrame:
    """
    TODO:
    - pentru fiecare medicament din graf calculați scorul de distanță
    - sortați crescător după distance
    """
    drugs = get_drug_nodes(B)
    results = []
    for d in drugs:
        score = compute_drug_disease_distance(B, d, disease_genes, mode)
        results.append({"drug": d, "distance": score})
    
    return pd.DataFrame(results).sort_values("distance").reset_index(drop=True)


# --------------------------
# Main
# --------------------------
if __name__ == "__main__":
    # TODO 1: verificați input-urile
    ensure_exists(DISEASE_GENES_TXT)

    # TODO 2: încărcați / construiți graful bipartit
    B_bipartite = load_bipartite_graph_or_build()

    # TODO 3: încărcați setul de disease genes
    disease_genes_set = load_disease_genes(DISEASE_GENES_TXT)

    # TODO 4: calculați ranking-ul medicamentelor după proximitate
    df_priority = rank_drugs_by_proximity(B_bipartite, disease_genes_set, mode="mean")

    # TODO 5: salvați rezultatele
    df_priority.to_csv(OUT_DRUG_PRIORITY, index=False)
    
    print("[INFO] Exercise 9.2 skeleton — completați TODO-urile și rulați scriptul.")