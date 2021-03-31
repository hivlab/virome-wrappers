__author__ = "Taavi Päll"
__copyright__ = "Copyright 2021, Taavi Päll"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"

import os

from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
prefix = os.path.splitext(os.path.basename(snakemake.output.fasta))[0]
assert (
    snakemake.output.qual == f"{prefix}.qual.txt"
), f"qual file name should be path/to/{prefix}.qual.txt"

shell(
    "(samtools mpileup -aa -A -d 0 -Q 0 {snakemake.input[0]} | ivar consensus -p {prefix} {extra}) {log}"
)
