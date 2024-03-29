__author__ = "Taavi Päll"
__copyright__ = "Copyright 2019, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

import pandas as pd
from Bio import SeqIO

if not isinstance(snakemake.input.coverage, str):
    raise ValueError("coverage must be file path")

if not isinstance(snakemake.input.contigs, str):
    raise ValueError("contigs must be file path")

# Filter contigs table by coverage
coverage = snakemake.params.get("avg_coverage", 2)

# Import coverage table
cov = pd.read_csv(snakemake.input.coverage, sep="\t")
good_coverage = cov["Avg_fold"] >= float(coverage)
cov_filtered = cov[good_coverage]
ids = list(cov_filtered["#ID"])
print("We have {} contigs".format(len(ids)))

# Subset contigs fasta
with open(snakemake.input.contigs, "rU") as input_handle, open(
    snakemake.output[0], "w"
) as output_handle:
    for contig in SeqIO.parse(input_handle, "fasta"):
        if contig.description in ids:
            SeqIO.write(contig, output_handle, "fasta")
