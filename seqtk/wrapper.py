__author__ = "Taavi Päll"
__copyright__ = "Copyright 2018, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

# Get parameters
seed = snakemake.params.get("seed", "11")
print("Running seqtk with seed:", seed)
frac = snakemake.params.get("frac", "1")
print("Sampling fraction of reads:", frac)

# Commands to run
if (frac > 0 and frac < 1):
    cmd = f"seqtk sample -s{seed} {{0}} {frac} > {{1}}"
    shell(cmd.format(snakemake.input[0], snakemake.output[0]))
    shell(cmd.format(snakemake.input[1], snakemake.output[1]))

else:
    cmd = "cp {0} {1}"
    shell(cmd.format(snakemake.input[0], snakemake.output[0]))
    shell(cmd.format(snakemake.input[1], snakemake.output[1]))
