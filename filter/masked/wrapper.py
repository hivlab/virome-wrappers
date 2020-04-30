__author__ = "Taavi Päll"
__copyright__ = "Copyright 2018, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from Bio import SeqIO

# Check if filtering parameters are supplied.
min_length = snakemake.params.get("min_length")
assert (
    min_length is not None
), "Please input a min_length parameter: length of consecutive sequence without N."
por_n = snakemake.params.get("por_n")
assert (
    por_n is not None
), "Please input a por_n parameter: % of total length of being masked."


def seq_filter(sequence, min_length, por_n):
    return (float(len(sequence)) - float(sequence.count("N"))) >= min_length and (
        float(sequence.count("N")) / float(len(sequence))
    ) * 100 <= por_n


def filter_N(masked, masked_filt, min_length, por_n, original=None, original_filt=None):
    """Filters out sequences with many N-s.
    Filters out sequeces with less than or equal to min_length non-N bases
    and sequences with more than por_n % N-s.
    """
    if original is not None and original_filt is not None:
        with open(masked_filt, "w") as masked_filt_handle, open(original_filt, "w") as original_filt_handle:
            for masked_record, original_record in zip(SeqIO.parse(masked, "fasta"), SeqIO.parse(original, "fasta")):
                if seq_filter(masked_record.seq, min_length, por_n):
                    SeqIO.write(masked_record, masked_filt_handle, "fasta")
                    SeqIO.write(original_record, original_filt_handle, "fasta")

    else:
        with open(masked_filt, "w") as masked_filt_handle:
            for record in SeqIO.parse(masked, "fasta"):
                if seq_filter(record.seq, min_length, por_n):
                    SeqIO.write(record, masked_filt_handle, "fasta")


if __name__ == "__main__":
    # Merge function arguments into dictionary.
    options = dict(snakemake.input)
    options.update(snakemake.output)
    options.update(snakemake.params)
    # Unwrap arguments and run filter.
    filter_N(**options)
