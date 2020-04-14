__author__ = "Ryan Dale"
__copyright__ = "Copyright 2016, Ryan Dale"
__email__ = "dalerr@niddk.nih.gov"
__license__ = "MIT"

import os
from snakemake.shell import shell
import tempfile

def file_name_prefix(file_path):
    file_name = os.path.basename(file_path)
    separator_index = file_name.index(".")
    return file_name[:separator_index]


_config = snakemake.params["fastq_screen_config"]
subset = snakemake.params.get("subset", 100000)
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell()

# snakemake.params.fastq_screen_config can be either a dict or a string. If
# string, interpret as a filename pointing to the fastq_screen config file.
# Otherwise, create a new tempfile out of the contents of the dict:
if isinstance(_config, dict):
    tmp = tempfile.NamedTemporaryFile(delete=False).name
    with open(tmp, "w") as fout:
        for label, indexes in _config["database"].items():
            print(indexes)
            fout.write(
                "\t".join(["DATABASE", label, indexes]) + "\n"
                )
    config_file = tmp
else:
    config_file = _config

# fastq_screen hard-codes filenames according to this prefix. We will send
# hard-coded output to a temp dir, and then move them later.
prefix = file_name_prefix(snakemake.input[0])
tempdir = tempfile.mkdtemp()

shell(
    "fastq_screen --outdir {tempdir} "
    "--force "
    "--aligner bwa "
    "--conf {config_file} "
    "--subset {subset} "
    "--threads {snakemake.threads} "
    "{extra} "
    "{snakemake.input[0]} "
    "{log}"
)

# Move output to the filenames specified by the rule
shell("mv {tempdir}/{prefix}_screen.txt {snakemake.output.txt}")
shell("mv {tempdir}/{prefix}_screen.png {snakemake.output.png}")

# Clean up temp
shell("rm -r {tempdir}")
if isinstance(_config, dict):
    shell("rm {tmp}")
