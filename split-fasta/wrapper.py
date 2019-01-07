__author__ = "Taavi Päll"
__copyright__ = "Copyright 2018, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from re import search
from math import ceil
from Bio import SeqIO

# http://biopython.org/wiki/Split_large_file
def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch

# Number of sequences in fasta file
fasta_file = open(snakemake.input[0])

seqs = sum([bool(search(r"^>", line)) for line in fasta_file])

# Calculate batch size given number of files
batch_size = ceil(seqs / snakemake.params[0])

# Split sequences into chunks based on batch size and write into files
record_iter = SeqIO.parse(snakemake.input[0], "fasta")

for n, batch in enumerate(batch_iterator(record_iter, batch_size), start = 1):
    output = [path for path in snakemake.output if search(r'[^0-9]{}[^0-9]+$'.format(n), path)]
    with open(output[0], "w") as handle:
        count = SeqIO.write(batch, handle, "fasta")
