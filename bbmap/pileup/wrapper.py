__author__ = "Taavi Päll"
__copyright__ = "Copyright 2019, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

def arg_c(args):
   """Concatenates input/output arguments with names."""
   argdict = dict(args)
   argdict.update((k, k + "=" + v) for k,v in argdict.items())
   return " ".join(list(argdict.values()))

# Get input/output and optional flags.
inputs = arg_c(snakemake.input)
outputs = arg_c(snakemake.output)
options = snakemake.params.get("options", "")

# Run command.
shell("pileup.sh"
      " {inputs}"
      " {outputs}"
      " {options}")
