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

# Setup log
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Temp files
temp_vcfgz = tempfile.TemporaryFile()
vcfgz_index = str(temp_vcfgz.name) + ".csi"
temp_consensus = tempfile.TemporaryFile()
temp_bed = tempfile.TemporaryFile()

shell(
    """
    (bgzip -c {snakemake.input.vcf} > {temp_vcfgz.name}
    bcftools index {temp_vcfgz.name}
    cat {snakemake.input.ref} | bcftools consensus {temp_vcfgz.name} > {temp_consensus.name}
    samtools view -h -b -q {mapping_quality} -f 0x3 {snakemake.input.alignment} | genomeCoverageBed -bga -ibam stdin | awk '$4 < {mask}' | bedtools merge > {temp_bed.name}
    bedtools maskfasta -fi {temp_consensus.name} -bed {temp_bed.name} -fo {snakemake.output.consensus}) {log}
    """
)

shell(
    """
    rm -f {temp_vcfgz.name}
    rm -f {vcfgz_index}
    rm -f {temp_consensus.name}
    rm -f {temp_bed.name}
    """
)
