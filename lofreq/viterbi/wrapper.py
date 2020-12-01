__author__ = "Taavi Päll"
__copyright__ = "Copyright 2020, Taavi Päll"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"


from snakemake.shell import shell
extra = snakemake.params.get("extra", "")

shell(
    "lofreq faidx {snakemake.input.ref} "
    "&& lofreq viterbi {extra} --ref {snakemake.input.ref} --out {snakemake.output[0]} {snakemake.input.bam}"
)
