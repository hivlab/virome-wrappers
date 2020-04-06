__author__ = "Taavi Päll"
__copyright__ = "Copyright 2021, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

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


inputs = snakemake.input


# Check for SE inputs and
n_files_diff = 0
if "in1" in inputs.keys():
    n_files = [len(v) for k, v in inputs.items() if k in ["in1", "in2"]]
    n_files_diff = n_files[0] - n_files[1]
    if n_files_diff < 0:
        raise Exception("Please check inputs: more elements in in2 than in in1.")


# Get inputs
inputs = arg_c(inputs, n_files_diff)

# Replace 'input' with 'in'
inputs = inputs.replace("input", "in")


# Pass only bbmap/bbwrap output parameters
clumpify_outputs = [
    "out",
    "out1",
    "out2",
]
outputs = dict(snakemake.output)
outputs = {k: v for k, v in outputs.items() if k in clumpify_outputs}

# Parse outputs to argument string
outputs = arg_c(outputs, 0)

# Get extra arguments
extra = snakemake.params.get("extra", "")

# Setup log
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    """
    (clumpify.sh {inputs} {outputs} reorder {extra}) {log}
    """
)
