__author__ = "Taavi Päll"
__copyright__ = "Copyright 2018, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

# Get parameters
seed = snakemake.params.get("seed", "11")
print("This is seed:", seed)
frac = snakemake.params.get("frac")
print("This is frac:", frac)
print(snakemake.params.frac)

# Command to run
if (frac > 0 and frac < 1):
    cmd = "seqtk sample -s{seed} {snakemake.input[%s]} {frac} > {snakemake.output[%s]}"
else:
    cmd = "ln -sr {snakemake.input[%s]} {snakemake.output[%s]}"

# Run command
for i in enumerate(snakemake.input):
    cmd = cmd % (i, i)
    shell(cmd)
