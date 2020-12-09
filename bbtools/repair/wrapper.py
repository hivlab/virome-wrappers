__author__ = "Taavi Päll"
__copyright__ = "Copyright 2020, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts


java_opts = get_java_opts(snakemake)


def parseIO(d):
    return " ".join([("in" if k == "input" else k) + "=" + v for k, v in d.items()])


inputs = parseIO(snakemake.input)
outputs = parseIO(snakemake.output)

# Get extra arguments
extra = snakemake.params.get("extra", "")

# Setup log
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    """
    (repair.sh {inputs} {outputs} {extra} {java_opts}) {log}
    """
)
