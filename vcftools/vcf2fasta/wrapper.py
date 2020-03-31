__author__ = "Taavi Päll"
__copyright__ = "Copyright 2020, Taavi Päll"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"


from snakemake.shell import shell


extra = snakemake.params.get("extra", "")


shell(
    "vcf2fasta {extra} -f {snakemake.input.ref} {snakemake.input.vcf} > {snakemake.output}"
)
