__author__ = "Taavi Päll"
__copyright__ = "Copyright 2018, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

import os
import pandas as pd
from Bio import SeqIO
from pandas.io.common import EmptyDataError

def touch(fname, times=None):
    """Emulates unix touch command.
    Source: https://stackoverflow.com/questions/1158076/implement-touch-using-python#1160227
    """
    with open(fname, 'a'):
        os.utime(fname, times)

def read_data(file):
    try:
        df = pd.read_table(file)
    except EmptyDataError:
        df = pd.DataFrame()
    return df

def parse_blast_result(blast_result, query, e_cutoff, outfmt, mapped, unmapped):
    """Finds out whether the BLAST best hit has a evalue lower than the cutoff. 
    If yes, outputs query information. If no, the sequence will be kept for further analysis. 
    Source: http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc93.
    Function expexts BLAST tabular format (outfmt 6).
    """
    # Import blast output table
    tab = read_data(blast_result)
    if len(tab.index) == 0:
        known_ids = set()
        touch(mapped)
    else:
        # Import column names, replace std when present, remove double quote
        std = "'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'"
        # Munge outfmt string
        outfmt = outfmt.replace("std", std)
        outfmt = outfmt.replace("'", "")
        # Get colnames from outfmt string
        colnames = [colname for colname in outfmt.split() if "6" not in colname]
        # Assign column names
        tab.columns = colnames
        # Filter results
        known = tab[(tab.evalue <= e_cutoff)]
        # Write seqs below threshold to file
        known.to_csv(mapped, sep = '\t', encoding = "utf-8", index = False)
        known_ids = set(known.qseqid)
    # Subset blast input
    with open(unmapped, "w") as out:
        for record in SeqIO.parse(str(query), "fasta"):
            if record.id not in known_ids:
                SeqIO.write(record, out, "fasta")

if __name__ == '__main__':
    options = dict(snakemake.input)
    options.update(snakemake.output)
    options.update(snakemake.params)
    # Unwrap arguments and run function
    parse_blast_result(**options)
