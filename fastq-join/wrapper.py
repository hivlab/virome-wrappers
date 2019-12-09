__author__ = "Taavi Päll"
__copyright__ = "Copyright 2018, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

# Check for correct number of input files
n = len(snakemake.input)
assert (
    n == 2 or n == 3
), "Input must contain 2 (paired-end) elements and optionally 'mate' input file."

# Get optional flags
options = snakemake.params.get("options", "")

# Setup log file
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Prepend -o flag output files and merge to a string
output = ["-o " + output for output in snakemake.output]
output = " ".join(output)

# Run command
shell("(fastq-join {options} {snakemake.input} {output}) {log}")
