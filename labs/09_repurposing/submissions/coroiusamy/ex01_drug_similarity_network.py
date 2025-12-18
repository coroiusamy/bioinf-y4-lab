"""
Exercise 9.1 — Drug–Gene Bipartite Network & Drug Similarity Network

Scop:
- să construiți o rețea bipartită drug–gene plecând de la un CSV
- să proiectați layer-ul de medicamente folosind similaritatea dintre seturile de gene
- să exportați un fișier cu muchiile de similaritate între medicamente

TODO:
- încărcați datele drug–gene
- construiți dict-ul drug -> set de gene țintă
- construiți graful bipartit drug–gene (NetworkX)
- calculați similaritatea dintre medicamente (ex. Jaccard)
- construiți graful drug similarity
- exportați tabelul cu muchii: drug1, drug2, weight
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Set, Tuple, List
import itertools
import networkx as nx
import pandas as pd
import pickle
import matplotlib.pyplot as plt

# --------------------------
# Config — adaptați pentru handle-ul vostru
# --------------------------
HANDLE = "coroiusamy"

# Input: fișier cu coloane cel puțin: drug, gene
DRUG_GENE_CSV = Path(f"data/work/{HANDLE}/lab09/drug_gene_{HANDLE}.csv")

# Output directory & files
OUT_DIR = Path(f"labs/09_repurposing/submissions/{HANDLE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_DRUG_SUMMARY = OUT_DIR / f"drug_summary_{HANDLE}.csv"
OUT_DRUG_SIMILARITY = OUT_DIR / f"drug_similarity_{HANDLE}.csv"
OUT_GRAPH_DRUG_GENE = OUT_DIR / f"network_drug_gene_{HANDLE}.pkl"

# Livrabile solicitate în PR
OUT_REPURPOSING_CSV = OUT_DIR / f"REPURPOSING_{HANDLE}.csv"
OUT_GRAPH_PNG = OUT_DIR / f"network_drug_gene_{HANDLE}.png"


def ensure_exists(path: Path) -> None:
    """
    TODO:
    - verificați că fișierul există
    - dacă nu, ridicați FileNotFoundError cu un mesaj clar
    """
    if not path.exists():
        raise FileNotFoundError(f"Fișierul de intrare nu a fost găsit: {path}")
    print(f"[INFO] Fișier găsit: {path}")


def load_drug_gene_table(path: Path) -> pd.DataFrame:
    """
    TODO:
    - citiți CSV-ul cu pandas
    - validați că există cel puțin coloanele: 'drug', 'gene'
    - returnați DataFrame-ul
    """
    df = pd.read_csv(path)
    if 'drug' not in df.columns or 'gene' not in df.columns:
        raise ValueError("DataFrame-ul trebuie să conțină coloanele 'drug' și 'gene'.")
    return df.dropna(subset=['drug', 'gene']).astype(str)


def build_drug2genes(df: pd.DataFrame) -> Dict[str, Set[str]]:
    """
    TODO:
    - construiți un dict: drug -> set de gene țintă
    - sugestie: folosiți groupby("drug") și aplicați set() pe coloana gene
    """
    return df.groupby("drug")["gene"].apply(set).to_dict()


def build_bipartite_graph(drug2genes: Dict[str, Set[str]]) -> nx.Graph:
    """
    TODO:
    - construiți graful bipartit:
      - nodurile 'drug' cu atribut bipartite="drug"
      - nodurile 'gene' cu atribut bipartite="gene"
      - muchii drug-gene
    """
    B = nx.Graph()
    for drug, gene_set in drug2genes.items():
        B.add_node(drug, bipartite="drug")
        for gene in gene_set:
            B.add_node(gene, bipartite="gene")
            B.add_edge(drug, gene)
    print(f"[INFO] Graf bipartit construit: {B.number_of_nodes()} noduri, {B.number_of_edges()} muchii.")
    return B


def summarize_drugs(drug2genes: Dict[str, Set[str]]) -> pd.DataFrame:
    """
    TODO:
    - construiți un DataFrame cu:
        drug, num_targets (numărul de gene țintă)
    - returnați DataFrame-ul
    """
    data = [{"drug": d, "num_targets": len(g)} for d, g in drug2genes.items()]
    df = pd.DataFrame(data)
    return df.sort_values(by="num_targets", ascending=False).reset_index(drop=True)


def jaccard_similarity(s1: Set[str], s2: Set[str]) -> float:
    """
    Calculați similaritatea Jaccard între două seturi de gene:
    J(A, B) = |A ∩ B| / |A ∪ B|
    """
    if not s1 and not s2: return 0.0
    intersection = len(s1 & s2)
    union = len(s1 | s2)
    return intersection / union if union > 0 else 0.0


def compute_drug_similarity_edges(
    drug2genes: Dict[str, Set[str]],
    min_sim: float = 0.0,
) -> List[Tuple[str, str, float]]:
    """
    TODO:
    - pentru toate perechile de medicamente (combinații de câte 2),
      calculați similaritatea Jaccard între seturile de gene
    - rețineți doar muchiile cu similaritate >= min_sim
    - returnați o listă de tuple (drug1, drug2, weight)
    """
    edges = []
    drugs = list(drug2genes.keys())
    for d1, d2 in itertools.combinations(drugs, 2):
        sim = jaccard_similarity(drug2genes[d1], drug2genes[d2])
        if sim >= min_sim:
            edges.append((d1, d2, sim))
    return edges


def edges_to_dataframe(edges: List[Tuple[str, str, float]]) -> pd.DataFrame:
    """
    TODO:
    - transformați lista de muchii (drug1, drug2, weight) într-un DataFrame
      cu coloanele: drug1, drug2, similarity
    """
    return pd.DataFrame(edges, columns=["drug1", "drug2", "similarity"])


# --------------------------
# Main
# --------------------------
if __name__ == "__main__":
    # TODO 1: verificați că fișierul de input există
    ensure_exists(DRUG_GENE_CSV)

    # TODO 2: încărcați tabelul drug-gene
    df_drug_gene = load_drug_gene_table(DRUG_GENE_CSV)

    # TODO 3: construiți mapping-ul drug -> set de gene
    drug2genes_map = build_drug2genes(df_drug_gene)

    # TODO 4: construiți graful bipartit și salvați-l (opțional)
    B_graph = build_bipartite_graph(drug2genes_map)
    with open(OUT_GRAPH_DRUG_GENE, 'wb') as f:
        pickle.dump(B_graph, f)
    
    # Salvare vizualizare PNG
    plt.figure(figsize=(12, 12))
    pos = nx.spring_layout(B_graph, k=0.1, iterations=15)
    nx.draw(B_graph, pos, node_size=5, alpha=0.3, with_labels=False)
    plt.title(f"Drug-Gene Network ({HANDLE})")
    plt.savefig(OUT_GRAPH_PNG)
    plt.close()

    # TODO 5: generați și salvați sumarul pe medicamente
    df_summary = summarize_drugs(drug2genes_map)
    df_summary.to_csv(OUT_DRUG_SUMMARY, index=False)
    # Livrabil 1
    df_summary.to_csv(OUT_REPURPOSING_CSV, index=False)

    # TODO 6: calculați similaritatea între medicamente
    sim_edges = compute_drug_similarity_edges(drug2genes_map, min_sim=0.1)
    df_similarity = edges_to_dataframe(sim_edges)
    df_similarity.to_csv(OUT_DRUG_SIMILARITY, index=False)

    print("[INFO] Exercise 9.1 skeleton — completați TODO-urile și rulați scriptul.")