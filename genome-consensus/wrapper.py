__author__ = "Taavi Päll"
__copyright__ = "Copyright 2020, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell
import tempfile
import os


# Parameters
mask = snakemake.params.get("mask", 5)
mapping_quality = snakemake.params.get("mapping_quality", 0)
name = snakemake.params.get("name", "\\1")


# Setup log
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Temp files
temp_vcfgz = tempfile.NamedTemporaryFile()
vcfgz_index = str(temp_vcfgz.name) + ".csi"
del_bed = tempfile.NamedTemporaryFile()
cov_bed = tempfile.NamedTemporaryFile()
mask_bed = tempfile.NamedTemporaryFile()
mask_bed = str(mask_bed.name) + ".bed"

shell(
    """
    (bgzip -c {snakemake.input.vcf} > {temp_vcfgz.name}
    bcftools index {temp_vcfgz.name}
    samtools view -h -b -q {mapping_quality} {snakemake.input.alignment} | genomeCoverageBed -bga -ibam stdin | awk '$4 < {mask}' | bedtools merge > {cov_bed.name}
    vcf2bed --deletions < {snakemake.input.vcf} > {del_bed.name}
    bedtools subtract -a {cov_bed.name} -b {del_bed.name} > {mask_bed}
    bcftools consensus -f {snakemake.input.ref} -m {mask_bed} {temp_vcfgz.name} | sed -E 's/^>(.*)/>{name}/' > {snakemake.output[0]}) {log}
    """
)
