__author__ = "Taavi Päll"
__copyright__ = "Copyright 2018, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

# Check inputs/arguments.
if not isinstance(snakemake.input, str) and len(snakemake.input) not in {1, 2}:
    raise ValueError("input must have 1 (single-end) or "
                     "2 (paired-end) elements")

# Extract arguments
options = snakemake.params.get("options", "")
html = snakemake.params.get("html", "fastp.html")
json = snakemake.params.get("json", "fastp.json")

# Setup log
log = snakemake.log_fmt_shell(stdout = False, stderr = True)

# Determine the number of input files
n = len(snakemake.input)
if n == 1:
    input = " -i {snakemake.input[0]} -o {snakemake.output[0]}"

elif n == 2: 
    input = (" -i {snakemake.input[0]} -o {snakemake.output[0]}"
             " -I {snakemake.input[1]} -O {snakemake.output[1]}")
else:
    raise ValueError("Unexpected number of inputs ({})".format(n))

# Run command
shell(
    "(fastp" + input + " {options}"
    " -h {html} -j {json}"
    " -w {snakemake.threads}) {log}")
