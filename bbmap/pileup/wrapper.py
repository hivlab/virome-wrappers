__author__ = "Taavi Päll"
__copyright__ = "Copyright 2019, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

# Function to concatenate arguments with names.
def arg_c(args):
   argdict = dict(args)
   arglist = ['{}={}'.format(k, ",".join(v if isinstance(v, type([])) else [v])) for k,v in argdict.items()]
   argstr = " ".join(arglist)
   return argstr

# Get arguments.
inputs = arg_c(snakemake.input)
outputs = arg_c(snakemake.output)
extra = snakemake.params.get("extra", "")

# Setup log.
log = snakemake.log_fmt_shell(stdout = False, stderr = True)

# Run command.
shell("(pileup.sh"
      " {inputs}"
      " {outputs}"
      " {options}) {log}")
