__author__ = "Taavi Päll"
__copyright__ = "Copyright 2019, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

# Check inputs/arguments.
inputs = dict(snakemake.input)
input_error_msg = "Input must contain following named elements: 'ref', 'in1', and 'in2'."
assert isinstance(inputs["ref"], str), "Must be path to reference fasta (e.g. contigs.fa)." 
input_names = list(inputs.keys())
assert input_names == ["ref", "in1", "in2"], input_error_msg

# Merge input paths with flags.
inputs.update((k, ",".join(v if isinstance(v, type([])) else [v])) for k,v in inputs.items())
input_flags = "ref={ref} in1 = {in1} in2 = {in2}".format(**inputs)

# Get optional flags.
options = snakemake.params.get("options", "")

# Run commands.
shell("bbwrap.sh"
      "{input_flags}"
      "out={snakemake.output.aln}"
      "{options}")
shell("pileup.sh"
      "in={snakemake.output.aln}"
      "out={snakemake.output.cov}")
