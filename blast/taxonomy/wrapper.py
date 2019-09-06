import pandas as pd
import subprocess
from ete3 import NCBITaxa
import numpy as np
import json
import sys
from itertools import islice

#chunk = int(sys.argv[1])
#print("We have chunk nr {}".format(chunk))
#"n_chunks = int(sys.argv[2])
#print("Total number of chunks is {}".format(n_chunks))
out_prefix="consensus_taxonomy"
n_chunks=3
chunk=1

# Create ncbi taxonomy query object
ncbi = NCBITaxa()

# Main taxonomy ranks
TAXONOMIC_RANKS = [ "no rank",
                    "species",
                    "genus",
                    "family",
                    "order",
                    "class",
                    "phylum",
                    "kingdom",
                    "superkingdom"]

# Dictionary that maps the rank to an integer,
# the lowest is "no rank": 0
NUMERIC_RANK = {TAXONOMIC_RANKS[i]: i for i in range(len(TAXONOMIC_RANKS))}

# Unknown rank
UNKNOWN_RANK = 'unknown'

# NCBI taxid representing unknown taxon
UNIDENTIFIED = 32644

blast = pd.read_csv("~/Downloads/queries.csv", index_col = "query", nrows = 50)
by_query = blast.groupby("query")
pp_sway = 1

def chunks(data, n_chunks = 100):
    it = iter(data)
    size = len(data) // n_chunks
    for i in range(0, len(data), size):
        yield {k:data[k] for k in islice(it, size)}

#def assign_consensus_taxonomy(blast_df, pp_sway, id, out_prefix):
consensus_taxonomy = []
for query, hits in by_query:
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
consensus_taxonomy.append({"query": query, "consensus": consensus[0], "pident": hits["pident"].aggregate("max"), "hits": hits.shape[0], "rank": json.dumps(ranks), "names": json.dumps(names)})

pd.DataFrame(consensus_taxonomy).to_csv("{}_{}.csv".format(out_prefix, chunk), index = False)

#for count, item in enumerate(chunks(by_query.indices, n_chunks), 1):
#    if count is chunk:
#        print(count, [*item])
#        print(blast.loc[[*item],])
#        assign_consensus_taxonomy(blast.loc[[*item],], pp_sway, count, "consensus_taxonomy")


