from Bio import SeqIO
import pandas as pd
    
def subset_records(virids, contigs, output):
    ids = pd.read_csv(virids).set_index("query", drop = False).index.to_list()
    with open(output, "w") as output_handle:
        for record in SeqIO.parse(contigs, "fasta"):
            if record.id.split()[0] in ids:
                SeqIO.write(record, output_handle, 'fasta')

subset_records(snakemake.input["virids"], snakemake.input["contigs"], snakemake.output[0])
