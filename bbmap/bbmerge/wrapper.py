__author__ = "Taavi Päll"
__copyright__ = "Copyright 2020, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell
import re

# Function to concatenate arguments with names
def arg_c(args, n):
    argdict = dict(args)
    # Convert values to list
    argdict = {k: (v if isinstance(v, type([])) else [v]) for k, v in argdict.items()}
    # Merge multiple argument values to comma separated str
    argdict = {
        k: (",".join(v + ["null"] * n_files_diff) if "in2" in k else ",".join(v))
        for k, v in argdict.items()
    }
    # Parse arguments for command line
    arglist = ["{}={}".format(k, v) for k, v in argdict.items()]
    argstr = " ".join(arglist)
    return argstr


assert len(snakemake.input) in [1, 2], "Input must have one or two files."

if len(snakemake.input) == 1:
    inputs = "in={}"
elif len(snakemake.input) == 2:
    inputs = "in1={} in2={}"

inputs = inputs.format(*snakemake.input)

# Pass only reformat outputs
reformat_outputs = [
    "out",
    "outu",
    "outu1",
    "outu2",
    "outinsert",
    "outadapter",
    "outc",
    "ihist",
]
outputs = dict(snakemake.output)
outputs = {k: v for k, v in outputs.items() if k in reformat_outputs}

# Parse outputs to argument string
outputs = arg_c(outputs, 0)

# Get extra arguments
extra = snakemake.params.get("extra", "")

# Setup log
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    """
    (bbmerge.sh {inputs} {outputs} {extra}) {log}
    """
)
