__author__ = "Taavi Päll"
__copyright__ = "Copyright 2020, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell
import re


def parseIO(d):
    return " ".join([("in" if k == "input" else k) + "=" + v for k, v in d.items()])


assert len(snakemake.input) in [1, 2], "Number of inputs must be one or two."


inputs = parseIO(snakemake.input)
outputs = parseIO(snakemake.output)

# Get extra arguments
extra = snakemake.params.get("extra", "")

# Setup log
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    """
    (filterbytile.sh {inputs} {outputs} {extra}) {log}
    """
)
