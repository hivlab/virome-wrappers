import pandas as pd
pd.options.mode.chained_assignment = None
from pkg_resources import parse_version
import warnings
from ete3 import NCBITaxa
import numpy as np
import argparse
import tarfile
import re

# Helper function to import tables
def safely_read_csv(path, **kwargs):
  try:
    return pd.read_csv(path, **kwargs)
  except pd.errors.EmptyDataError:
    pass

class pandasVersionWarning(UserWarning):
    pass

if parse_version(pd.__version__) < parse_version("0.24.0"):
    warnings.warn("Pandas version is {} and taxon ids will be floating point numbers in results\nand column order is not preserved. To avoid this, please install pandas>=0.24.0.".format(pd.__version__), pandasVersionWarning)

# Create BLAST BD and query classes
class BlastDB:

    def __init__(self, dbfile=None):
        self.ncbi = NCBITaxa(dbfile=dbfile)
    
    def get_lineage(self, taxid):
        return self.ncbi.get_lineage(taxid)

    def get_rank(self, taxids):
        return self.ncbi.get_rank(taxids)
    
    def get_topology(self, lineage):
        return self.ncbi.get_topology(lineage)
        
    def translate_to_names(self, taxids):
        return self.ncbi.translate_to_names(taxids)

class BlastTaxonomy(BlastDB):
    
    def __init__(self, results, query_key = "qseqid", taxid_key = "staxid", pp_sway=1, ranks_of_interest=None, dbfile=None):
        
        BlastDB.__init__(self, dbfile)
        self.query_key = query_key
        self.taxid_key = taxid_key
        self.by_query = results.groupby(self.query_key)
        self.pp_sway = pp_sway
        if ranks_of_interest:
            self.ranks_of_interest = ranks_of_interest
        else:
            self.ranks_of_interest = ["superkingdom", "order", "family", "genus", "species"]
        self.unidentified = 32644
    
    def get_consensus_taxonomy(self):
        consensus_taxonomy = []
        for query, hits in self.by_query:
            if hits.shape[0] > 1:
                # Keep only one top hit from each taxon
                hits["pident_rank"] = hits.groupby([self.taxid_key])["pident"].rank(method = "first", ascending = False)
                hits = hits[hits["pident_rank"] == 1]
                # Try to remove unidentified taxa
                hits["name"] = hits[self.taxid_key].apply(lambda x: self.translate_to_names([x])[0])
                unidentified = hits["name"].apply(lambda x: bool(re.search("unident", x)))
                identified = np.invert(unidentified)
                if sum(identified) >= 1:
                    # Keeping only identified  taxids
                    hits = hits[identified]
                    # Filtering by percent identity
                    pident_threshold = hits["pident"].aggregate("max") - self.pp_sway
                    within = hits["pident"].apply(lambda x: x >= pident_threshold)
                    hits_filtered = hits[within]
                    # Getting consensus taxonomy
                    taxlist = hits_filtered[self.taxid_key].tolist()
                    if len(taxlist) > 1:
                        tree = self.get_topology(taxlist)
                        for node in tree.traverse():
                            if node.is_root():
                                consensus = [node.name]
                    else:
                        consensus = taxlist
                else:
                    consensus = [self.unidentified]
            else:
                consensus = hits[self.taxid_key].tolist()   
            lineage = self.get_lineage(*consensus)
            ranks = self.get_rank(lineage)
            con_lin = {rank: id for id, rank in ranks.items() if rank in self.ranks_of_interest}
            lin_names = {rank + "_name": self.translate_to_names([id])[0] for rank, id in con_lin.items()}
            consensus_taxonomy.append(dict({"query": query, "consensus": consensus[0], "pident": hits["pident"].aggregate("max"), "hits": hits.shape[0]}, **con_lin, **lin_names))
        consensus_taxonomy = pd.DataFrame(consensus_taxonomy)
        # Convert tax_ids to integers
        consensus_taxonomy[self.ranks_of_interest] = consensus_taxonomy[self.ranks_of_interest].apply(lambda x: pd.Series(x, dtype = "Int64"))
        return consensus_taxonomy

def blast_taxonomy(input, output, sep = "\t", **kwargs):

    # Import file with BLAST results
    run = []
    for file in input:
        if tarfile.is_tarfile(file):
            with tarfile.open(file, "r:*") as tar:
                splits = []
                for member in tar.getmembers():
                    m = tar.extractfile(member)
                    splits.append(safely_read_csv(m, sep = sep))
                run.append(pd.concat(splits))
        else:
            run.append(safely_read_csv(file, sep = sep))

    if all(v is None for v in run):
        consensus_taxonomy = pd.DataFrame()
    else:
        results = pd.concat(run, sort = False)
        bt = BlastTaxonomy(results, **kwargs)
        consensus_taxonomy = bt.get_consensus_taxonomy()
    
    with open(output, "w") as outfile:
        consensus_taxonomy.to_csv(outfile, index = False)

if __name__ == "__main__":

    blast_taxonomy(snakemake.input, snakemake.output[0], **snakemake.params)

    
