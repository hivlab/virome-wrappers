__author__ = "Taavi Päll"
__copyright__ = "Copyright 2020, Taavi Päll"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"


from snakemake.shell import shell


# Parse params
extra = snakemake.params.get("extra", "")


shell("samtools rmdup {extra} {snakemake.input[0]} {snakemake.output[0]}")