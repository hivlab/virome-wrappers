__author__ = "Taavi Päll"
__copyright__ = "Copyright 2018, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from Bio import SeqIO

# Check if filtering parameters are supplied.
min_length = snakemake.params.get("min_length")
assert min_length is not None, "Please input a min_length parameter: length of consecutive sequence without N."
por_n = snakemake.params.get("por_n")
assert por_n is not None, "Please input a por_n parameter: % of total length of being masked."

def filter_N(masked, masked_filt, min_length, por_n, original = None, original_filt = None):
  """Filters out sequences with many N-s.
  Filters out sequeces with less than or equal to min_length non-N bases
  and sequences with more than por_n % N-s.
  """
  ma = SeqIO.parse(str(masked), "fasta")
  if (original is not None and original_filt is not None):
    original = SeqIO.index(str(original), "fasta")
    with open(masked_filt, "w") as maf, open(original_filt, "w") as orf:
      for record in ma:
        sequence = str(record.seq).upper()
        if ((float(len(sequence)) - float(sequence.count("N"))) >= min_length and (float(sequence.count("N")) / float(len(sequence))) * 100 <= por_n):
          SeqIO.write(record, maf, "fasta")
          SeqIO.write(original[record.id], orf, "fasta")
  else:
    with open(masked_filt, "w") as maf:
      for record in ma:
        sequence = str(record.seq).upper()
        if ((float(len(sequence)) - float(sequence.count("N"))) >= min_length and (float(sequence.count("N")) / float(len(sequence))) * 100 <= por_n):
          SeqIO.write(record, maf, "fasta")

if __name__ == '__main__':
  # Merge function arguments into dictionary.
  options = dict(snakemake.input)
  options.update(snakemake.output)
  options.update(snakemake.params)
  # Unwrap arguments and run filter.
  filter_N(**options)
