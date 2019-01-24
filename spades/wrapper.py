__author__ = "Taavi Päll"
__copyright__ = "Copyright 2019, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell
from os.path import dirname

# Check inputs/arguments.
if not isinstance(snakemake.input, str) and len(snakemake.input) not 2:
    raise ValueError("input must have 2 (paired-end) elements")

# Extract arguments
options = snakemake.params.get("options", "")

# Get output dir name from output path where spades writes its output files
# See spades manual http://cab.spbu.ru/files/release3.12.0/manual.html for SPAdes output
# Parse output to a list just in case there is more than one output file specified,
# pick output dir from the first output file path
output = list(snakemake.output.values())
output_dir = dirname(output[0])

# Setup log
log = snakemake.log_fmt_shell(stdout = False, stderr = True)

shell("mkdir -p {output_dir}")
shell("spades.py {options} -1 {snakemake.input[0]} -2 {snakemake.input[1]} -o {output_dir} 2> {log}")
