__author__ = "Taavi Päll"
__copyright__ = "Copyright 2019, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

# Get parameters
extra = snakemake.params.get("extra", "")

# Run tantan
shell("tantan {extra} {snakemake.input} > {snakemake.output}")
