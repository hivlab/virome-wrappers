__author__ = "Taavi Päll"
__copyright__ = "Copyright 2019, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

# Parse inputs/arguments.
def arg_c(args):
   """Concatenates input/output arguments with names."""
   argdict = dict(args)
   argdict.update((k, k + "=" + v if len(k) > 0 else "in=" + v) for k,v in argdict.items())
   return " ".join(list(argdict.values()))

# Get input/output and optional flags.
if len(snakemake.input) == 1:
   inputs = "in=" + str(snakemake.input)
else:
   inputs = arg_c(snakemake.input)


print("snakemake.input is " + str(snakemake.input))
print("Length of snakemake.input is " + len(snakemake.input))
print("Print out inputs:")
print(inputs)

print("Asserting...")
assert len(inputs) > 0, "Input error. Input can have only one unnamed value assigned to variable 'in'. All other inputs must be named."
outputs = arg_c(snakemake.output)
options = snakemake.params.get("options", "")

# Setup log.
log = snakemake.log_fmt_shell(stdout = False, stderr = True)

# Run command.
shell("(pileup.sh"
      " {inputs}"
      " {outputs}"
      " {options}) {log}")
