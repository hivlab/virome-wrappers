__author__ = "Taavi Päll"
__copyright__ = "Copyright 2020, Taavi Päll"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"


from snakemake.shell import shell


indexfeature = snakemake.params.get("indexfeature", "")
refmaker = snakemake.params.get("refmaker", "")


shell(
        "gatk IndexFeatureFile {indexfeature} -I {snakemake.input.vcf} -O {snakemake.output.idx}"
    )

shell(
        "gatk FastaAlternateReferenceMaker {refmaker} -V {snakemake.input.vcf} -R {snakemake.input.ref} -O {snakemake.output.fasta}"
    )
