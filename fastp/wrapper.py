__author__ = "Taavi Päll"
__copyright__ = "Copyright 2018, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

# Extract arguments
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout = False, stderr = True)

# Run command
shell(
    "(fastp -i {snakemake.input[0]} -I {snakemake.input[1]}"
    " -o {snakemake.output.pair1} -O {snakemake.output.pair2}"
    " {extra}"
    " -h {snakemake.output.html} -j {snakemake.output.json}"
    " -w {snakemake.threads})"
    " {log}")
