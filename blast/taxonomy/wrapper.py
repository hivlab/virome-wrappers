import pandas as pd
import subprocess
from ete3 import NCBITaxa
import numpy as np
import json
import sys
from itertools import islice
import argparse

parser = argparse.ArgumentParser(description = "Process BLAST+ taxonomy in chunks.")
parser.add_argument("--infile", dest = "filename", type = argparse.FileType("r"), 
                required=True,
                   help = "Path to input file with BLAST results")
parser.add_argument("-i", dest = "index", type = int, nargs = 1,
                    required=True,
                   help = "Array index, an integer")
parser.add_argument("-s", dest = "size", type = int, nargs = 1,
                   required = True, help = "Array size, an integer")
parser.add_argument("--nrows", dest = "nrows", default = None,
                    help = "Number of rows of file to read. Useful for reading pieces of large files.")
args = parser.parse_args()

# Output file prefix
out_prefix = "consensus_taxonomy"

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

# Unknown rank and unidentified taxid
UNKNOWN_RANK = 'unknown'
UNIDENTIFIED = 32644

# Percent point of difference from maximum identity hit
pp_sway = 1

# Taxonomic ranks of interest
RANKS_OF_INTEREST = ['superkingdom', 'order', 'family', 'genus', 'species']

def split(x, n):
    k, m = divmod(len(x), n)
    return (x[i * k + min(i, m): (i + 1) * k + min(i + 1, m)] for i in range(n))

def get_normalised_lineage(taxid, ranks_of_interest, taxonomic_ranks):
    print(taxid)
    # Dictionary that maps taxonomic ranks to integers,
    # the lowest is "no rank": 0
    numeric_rank = {taxonomic_ranks[i]: i for i in range(len(taxonomic_ranks))}
    unidentified = 32644

    # Get lineage for taxid and map it to numeric rank
    lineage = ncbi.get_lineage(taxid)
    _, rank_of_query = ncbi.get_rank([lineage[-1]]).popitem()
    num_rank_of_query = numeric_rank[rank_of_query]
    ranks = ncbi.get_rank(lineage)
    invert_dict = {v: k for k, v in ranks.items()}
    lineage_ranks = invert_dict.keys()
    
    # Get ranks and taxids
    # if rank is not present in lineage, say is unidentified
    normalised_lineage = dict()
    for rank in ranks_of_interest:
        num_rank_to_get = numeric_rank[rank]
        # must be above taxid in hierarchy
        if num_rank_to_get >= num_rank_of_query:
            if rank in lineage_ranks:
                normalised_lineage[rank] = invert_dict[rank]
            else:
                normalised_lineage[rank] = unidentified

    return(normalised_lineage)

def assign_consensus_taxonomy(blast_df, pp_sway, out_prefix, **kwargs):
    consensus_taxonomy = []
    for query, hits in blast_df.groupby("query"):
        if hits.shape[0] > 1:
            hits = hits.sort_values(by = "evalue")
            pident_threshold = hits["pident"].aggregate("max") - pp_sway
            within = hits["pident"].apply(lambda x: x >= pident_threshold)
            hits_filtered = hits[within]
            taxlist = hits_filtered["tax_id"].tolist()
            if len(taxlist) > 1:
                lineage = []
                for tax in taxlist:
                    normalised_lineage = get_normalised_lineage(tax, **kwargs)
                    rev_normalised_lineage = {v: k for k, v in normalised_lineage.items()}
                    lineage.append(set(rev_normalised_lineage))
                lineage_intersect = list(set.intersection(*lineage))
                root_tree = ncbi.get_topology(lineage_intersect)
                consensus = root_tree.get_leaf_names()
            else:
                consensus = taxlist
        else:
            consensus = hits["tax_id"].tolist()
        con_lin = get_normalised_lineage(consensus[0], **kwargs)
        consensus_taxonomy.append(dict({"query": query, "consensus": consensus[0], "pident": hits["pident"].aggregate("max"), "hits": hits.shape[0]}, **con_lin))
    pd.DataFrame(consensus_taxonomy).to_csv("{}_{}.csv".format(out_prefix, args.index[0]), index = False)

# Import parsed and evalue filtered BLAST+ results
blast = pd.read_csv(args.filename, index_col = "query", nrows = int(args.nrows) if args.nrows else args.nrows)
by_query = blast.groupby("query")
queries = [*by_query.indices]

# Split queries into n chunks
query_chunks = list(split(queries, args.size[0]))

# Process one chunk
queries_to_process = query_chunks[args.index[0] - 1]
assign_consensus_taxonomy(blast.loc[queries_to_process], pp_sway, out_prefix, ranks_of_interest = RANKS_OF_INTEREST, taxonomic_ranks = TAXONOMIC_RANKS)
