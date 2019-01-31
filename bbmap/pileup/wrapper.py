__author__ = "Taavi Päll"
__copyright__ = "Copyright 2019, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

# Check inputs/arguments.
def arg_c(args):
   """Concatenates input/output arguments with names."""
   argdict = dict(args)
   argdict.update((k, k + "=" + v) for k,v in argdict.items())
   return " ".join(list(argdict.values()))

# Get input/output and optional flags.
# First see if we have unnamed values to 'in' argument.
args = dict(snakemake.input)
unnamed = [v for v in snakemake.input if v not in set(args.values())]

# Get named arguments.
inputs = arg_c(snakemake.input)

# Parse value for 'in' argument and merge with other named arguments. 
assert len(unnamed) <= 1, "More than one value to 'in' argument."
if len(unnamed) == 1:
   inputs_to_in = "in={}".format(",".join(unnamed))
   inputs = " ".join([inputs_to_in, inputs])

outputs = arg_c(snakemake.output)
options = snakemake.params.get("options", "")

# Setup log.
log = snakemake.log_fmt_shell(stdout = False, stderr = True)

# Run command.
shell("(pileup.sh"
      " {inputs}"
      " {outputs}"
      " {options}) {log}")
