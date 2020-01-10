__author__ = "Taavi Päll"
__copyright__ = "Copyright 2019, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell
import re

# bbduk adapter trimming parameters.
bbduk = snakemake.params.get("bbduk", "")

# Check and fix if int=f flag is not set in bbduk params
if "int=f" not in bbduk:
    bbduk = bbduk + " int=f"

# Subsampling parameters.
frac = snakemake.params.get("frac", 1.0)
seed = snakemake.params.get("seed", 11)
print("Sampling {} fraction of reads & using seed {}.".format(frac, seed))

# Parse stdin
if re.search("fastq$|fq$", snakemake.output.merged):
    stdin = "stdin.fq"
elif re.search("fastq.gz$|fq.gz$", snakemake.output.merged):
    stdin = "stdin.fq.gz"
else:
    raise ValueError("Merge output must be fastq (.fastq/.fq).")

# Preprocessing command to run.
commands = [
    "bbmerge.sh in1={snakemake.input[0]} in2={snakemake.input[1]} outa={snakemake.output.adapters}",
    "bbmerge.sh in1={snakemake.input[0]} in2={snakemake.input[1]} out={snakemake.output.merged} outu={snakemake.output.unmerged} adapters={snakemake.output.adapters}",
    "cat {snakemake.output.merged} {snakemake.output.unmerged} | bbduk.sh in={stdin} out={snakemake.output.trimmed} ref={snakemake.output.adapters} {bbduk}",
]

# Run preprocessing commands.
for cmd in commands:
    shell(cmd)

# If sample fraction is given, subsample reads using seed.
# Otherwise copy trimmed reads to final output.
if frac and float(frac) < 1.0:
    shell(
        "reformat.sh in={snakemake.output.trimmed} out={snakemake.output.sampled} samplerate={frac} sampleseed={seed}"
    )
else:
    shell("cp {snakemake.output.trimmed} {snakemake.output.sampled}")
