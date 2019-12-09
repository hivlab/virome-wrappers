from Bio import SeqIO

# https://www.biostars.org/p/10162/
def get_ids(source):
    ids = map(lambda x: x.id, source)
    return set(ids)


def subset_records(source, ids):
    for record in source:
        if record.id in ids:
            yield record


# Parse fasta
fasta = SeqIO.parse(snakemake.input[0], "fasta")

# Get fasta ids
fasta_ids = get_ids(fasta)

# Subset another fasta using fasta_ids
subset = subset_records(SeqIO.parse(snakemake.input[1], "fasta"), fasta_ids)

# Write subset to a file
subset_count = SeqIO.write(subset, snakemake.output[0], "fasta")
