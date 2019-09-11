__author__ = "Taavi Päll"
__copyright__ = "Copyright 2019, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

# Setup log
log = snakemake.log_fmt_shell(stdout = False, stderr = True)

# Merge function arguments into dictionary.
# Look for available options for the BLAST+ command-line applications
# https://www.ncbi.nlm.nih.gov/books/NBK279684/
options = dict(snakemake.input)
options.update(snakemake.output)
options.update(snakemake.params)

# Remove BLAST+ program from dictionary before options formatting
program = options.pop("program")

options = " ".join(["-{} {}".format(k, v) for k, v in options.items()])

shell("({program}"
      " {options}"
      " -num_threads {snakemake.threads}) {log}")
