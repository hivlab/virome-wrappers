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

# Get inputs.
inputs = arg_c(snakemake.input)

# Replace 'input' with 'in'
inputs = inputs.replace("input", "in")

# Get outputs.
outputs = arg_c(snakemake.output)

# Pass only bbmap/bbwrap output parameters
bbmap_outputs = ["out", "outu", "outm", "scafstats", "refstats", "bhist", "qhist", "aqhist", "bqhist", "lhist", "ihist", "ehist", "qahist", "indelhist", "mhist", "gchist", "idhist", "covstats", "rpkm", "covhist", "basecov", "bincov"]
outputs = {k:v for k, v in outputs.items() if k in bbmap_outputs}

# Get extra arguments.
extra = snakemake.params.get("extra", "")

# Setup log.
log = snakemake.log_fmt_shell(stdout = False, stderr = True)

# Run commands.
shell("(bbwrap.sh"
      " {inputs}"
      " {outputs}"
      " {extra}) {log}")
