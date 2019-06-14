from Bio import SeqIO
import pandas as pd
    
def subset_records(source, ids):
    for record in source:
        if record.id.split()[0] in ids:
            yield record
    
virids = pd.read_csv(snakemake.input.virids[0]).set_index("query", drop = False).index.to_list()
subset = subset_records(SeqIO.parse(snakemake.input.contigs[0], "fasta"), virids)
with open(snakemake.output, "w") as out:
    SeqIO.write(subset, out, 'fasta')
