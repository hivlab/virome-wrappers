__author__ = "Taavi Päll"
__copyright__ = "Copyright 2019, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

inputs = dict(snakemake.input)
outputs = dict(snakemake.output)
options = snakemake.params.get("options", "")

# Merge input and output paths with flags.
inputs.update((k, k + "=" + v) for k,v in inputs.items())
input_files = " ".join(list(inputs.values()))
outputs.update((k, k + "=" + v) for k,v in outputs.items())
output_files = " ".join(list(inputs.values()))

# Run command.
shell("pileup.sh"
      " {input_files}"
      " {output_files}"
      " {options}")
