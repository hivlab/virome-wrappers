__author__ = "Taavi Päll"
__copyright__ = "Copyright 2020, Taavi Päll"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"


from snakemake.shell import shell
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    "lofreq alnqual {extra}"
    " {snakemake.input.bam}"
    " {snakemake.input.ref}"
    " > {snakemake.output[0]} {log}"
)
