import pandas as pd
import subprocess
from ete3 import NCBITaxa
from sklearn.cluster import SpectralClustering
import numpy as np
import json

ncbi = NCBITaxa()

blast = pd.read_csv("~/Downloads/queries.csv", index_col = "query", nrows = 100)
print(blast)
by_query = blast.groupby("query")
pp_sway = 1

cons_tax = []
with open("consensus_taxonomy.csv", )
for query, hits in by_query:
    print(query)
    if hits.shape[0] > 1:
        hits = hits.sort_values(by = "evalue")
        pident_threshold = hits["pident"].aggregate("max") - pp_sway
        within = hits["pident"].apply(lambda x: x >= pident_threshold)
        hits_filtered = hits[within]
        taxlist = hits_filtered["tax_id"].tolist()
        lineage = []
        for tax in taxlist:
            lineage.append(set(ncbi.get_lineage(tax)))
            lineage_intersect = list(set.intersection(*lineage))
        root_tree = ncbi.get_topology(lineage_intersect)
        consensus = root_tree.get_leaf_names()
        assert(len(consensus) == 1)
    else:
        consensus = hits["tax_id"].tolist()
    con_lin = []
    for tax in consensus:
        con_lin.append(ncbi.get_lineage(tax))
    ranks = []
    names = []
    for tax in con_lin:
        ranks.append(ncbi.get_rank(tax))
        names.append(ncbi.get_taxid_translator(tax))
    df = pd.DataFrame([{"query": query, "consensus": consensus[0], "hits": hits.shape[0], "rank": json.dumps(ranks), "names": json.dumps(names)}])
    df.to_csv()

