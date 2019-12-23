__author__ = "Taavi Päll"
__copyright__ = "Copyright 2019, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell


# Function to concatenate arguments with names
def arg_c(args):
    argdict = dict(args)
    arglist = [
        "{}={}".format(k, ",".join(v if isinstance(v, type([])) else [v]))
        for k, v in argdict.items()
    ]
    argstr = " ".join(arglist)
    return argstr


# Get inputs
inputs = arg_c(snakemake.input)

# Replace 'input' with 'in'
inputs = inputs.replace("input", "in")

# Check for SE inputs and
if "in1" in inputs.keys():
    n_files = [len(v) for k, v in inputs.items() if k in ["in1", "in2"]]
    n_files_diff = n_files[0] - n_files[1]
    if n_files_diff < 0:
        raise Exception("Please check inputs: more elements in in2 than in in1.")
    else:
        inputs = {
            k: (",".join(v + ["null" * n_files_diff]) if "in2" in k else ",".join(v))
            for k, v in inputs.items()
        }


# Pass only bbmap/bbwrap output parameters
bbmap_outputs = [
    "out",
    "outu",
    "outm",
    "scafstats",
    "refstats",
    "bhist",
    "qhist",
    "aqhist",
    "bqhist",
    "lhist",
    "ihist",
    "ehist",
    "qahist",
    "indelhist",
    "mhist",
    "gchist",
    "idhist",
    "covstats",
    "rpkm",
    "covhist",
    "basecov",
    "bincov",
]
outputs = dict(snakemake.output)
outputs = {k: v for k, v in outputs.items() if k in bbmap_outputs}

# Parse outputs to argument string
outputs = arg_c(outputs)

# Get extra arguments
extra = snakemake.params.get("extra", "")

# Setup log
log = snakemake.log_fmt_shell(stdout=False, stderr=True)


# Run commands
shell("(bbwrap.sh {inputs} {outputs} {extra}) {log}")
if "bamscript=bs.sh" in extra:
    shell("./bs.sh")
