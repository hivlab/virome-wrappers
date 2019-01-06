__author__ = "Taavi Päll"
__copyright__ = "Copyright 2018, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

seed = snakemake.params.get("seed", "11")
frac = snakemake.params.get("frac", "1")

if (frac > 0 and frac < 1):
    shell("seqtk sample -s{seed} {snakemake.input[0]} {frac} > {snakemake.output[0]}")
    shell("seqtk sample -s{seed} {snakemake.input[1]} {frac} > {snakemake.output[1]}")
else:
    shell("ln -sr {snakemake.input[0]} {snakemake.output[0]}")
    shell("ln -sr {snakemake.input[1]} {snakemake.output[1]}")
