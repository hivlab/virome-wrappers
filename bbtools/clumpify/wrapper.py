__author__ = "Taavi Päll"
__copyright__ = "Copyright 2020, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell


def parseIO(d):
    return " ".join([("in" if k == "input" else k) + "=" + v for k, v in d.items()])


assert len(snakemake.input) == len(
    snakemake.output
), "Number of inputs and outputs must be equal (one or two)."


inputs = parseIO(snakemake.input)
outputs = parseIO(snakemake.output)

# Get extra arguments
extra = snakemake.params.get("extra", "")

# Setup log
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    """
    (clumpify.sh {inputs} {outputs} reorder {extra}) {log}
    """
)
