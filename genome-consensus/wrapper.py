__author__ = "Taavi Päll"
__copyright__ = "Copyright 2020, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell


# Get extra arguments
mask = snakemake.params.get("mask", 20)
extra = snakemake.params.get("extra", "")
reads=snakemake.input.reads
reads=",".join([reads] if isinstance(reads, str) else reads)

# Setup log
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    """
    (bgzip -c {snakemake.input.vcf} > {snakemake.output.vcfgz}
    bcftools index {snakemake.output.vcfgz}
    cat {snakemake.input.ref} | bcftools consensus {snakemake.output.vcfgz} > {snakemake.output.consensus}
    bbwrap.sh ref={snakemake.output.consensus} in={reads} out={snakemake.output.sam} nodisk append {extra}
    samtools view -h -b -q 20 -f 0x3 {snakemake.output.sam} | genomeCoverageBed -bga -ibam stdin | awk '$4 < {mask}' | bedtools merge > {snakemake.output.bed}
    bedtools maskfasta -fi {snakemake.output.consensus} -bed {snakemake.output.bed} -fo {snakemake.output.consensus_masked}) {log}
    """
)
