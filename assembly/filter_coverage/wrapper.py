__author__ = "Taavi Päll"
__copyright__ = "Copyright 2019, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

import pandas as pd
from Bio import SeqIO

# Filter contigs table by coverage
cov = pd.read_table(snakemake.input.coverage)
good_coverage = cov["Avg_fold"] >= snakemake.params.avg_coverage
cov_filtered = cov[good_coverage]
ids = list(cov_filtered["#ID"])

# Subset contigs fasta
with open(snakemake.input.contigs, 'rU') as input_fasta, open(snakemake.output, 'w') as filtered_fasta:
    for contig in SeqIO.parse(input_fasta, 'fasta'):
        if contig.name in ids:
            SeqIO.write(contig, filtered_fasta, 'fasta')
