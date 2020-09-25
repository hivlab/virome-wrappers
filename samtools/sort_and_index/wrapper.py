__author__ = "Taavi Päll"
__copyright__ = "Copyright 2020, Taavi Päll"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"


import os
from snakemake.shell import shell


prefix = os.path.splitext(snakemake.output.sorted)[0]


threads = "" if snakemake.threads <= 1 else f" -@ {snakemake.threads - 1} "


shell(
    """
    samtools sort {snakemake.params} {threads} -o {snakemake.output.sorted} -T {prefix} {snakemake.input[0]} \
    && samtools index {snakemake.output.sorted} {snakemake.output.index}
    """
)
