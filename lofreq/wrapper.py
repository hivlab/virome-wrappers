__author__ = "Taavi Päll"
__copyright__ = "Copyright 2020, Taavi Päll"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"


from snakemake.shell import shell
extra = snakemake.params.get("extra", "")

shell(
    "lofreq call-parallel {extra} --pp-threads {snakemake.threads} -f {snakemake.input.ref} -o {snakemake.output[0]} {snakemake.input.bam}"
)
