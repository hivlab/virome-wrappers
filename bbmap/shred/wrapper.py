__author__ = "Taavi Päll"
__copyright__ = "Copyright 2019, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

# Get extra arguments
extra = snakemake.params.get("extra", "")

# Setup log
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Run commands
shell("(shred.sh in={snakemake.input[0]} out={snakemake.output[0]} {extra}) {log}")
