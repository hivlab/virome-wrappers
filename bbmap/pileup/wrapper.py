__author__ = "Taavi Päll"
__copyright__ = "Copyright 2019, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

# Check inputs/arguments.
def arg_c(args):
   """Concatenates input/output arguments with names."""
   argdict = dict(args)
   argdict.update((k, k + "=" + v if len(k) > 0 else "in=" + v) for k,v in argdict.items())
   return " ".join(list(argdict.values()))

# Check that input has only max one unnamed argument.
argdict = dict(snakemake.input)
assert sum([len(k) == 0 for k in list(argdict.keys())]) <= 1, "Unnamed input is reserved for 'in' argument. Please see pileup.sh help for available arguments."

# Get input/output and optional flags.
inputs = arg_c(snakemake.input)
print(inputs)
outputs = arg_c(snakemake.output)
print(outputs)
options = snakemake.params.get("options", "")

# Setup log.
log = snakemake.log_fmt_shell(stdout = False, stderr = True)

# Run command.
shell("(pileup.sh"
      " {inputs}"
      " {outputs}"
      " {options}) {log}")
