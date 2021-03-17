__author__ = "Taavi Päll"
__copyright__ = "Copyright 2021, Taavi Päll"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"

from snakemake.shell import shell


def arg_c(args):
    return " ".join([f"--{k} {v}" for k, v in args.items()])


# Parse parameters
task = snakemake.params.get("task", None)
assert (
    task is not None
), "Task is a required parameter. Must be one of 'uchime_denovo', 'uchime2_denovo', 'uchime3_denovo' or 'uchime_ref'."
extra = snakemake.params.get("extra", "")

# Parse output
output = snakemake.output
output = {
    k: v
    for k, v in output.items()
    if k in ["chimeras", "nonchimeras", "uchimealns", "uchimeout"]
}


# Parse input
assert len(snakemake.input) in [
    1,
    2,
], "Input must be one fasta file OR fastafile=fastafile and db=fastafile."
if len(snakemake.input) == 1:
    input = snakemake.input[0]
else:
    assert all(
        [k in snakemake.input.keys() for k in ["fastafile", "db"]]
    ), "Please provide two named input arguments: fastafile=fastafile & db=fastafile"
    input = f"{snakemake.input['fastafile']} --db {snakemake.input['db']}"


# Parse outputs to argument string
output = arg_c(output)

# Get extra arguments
extra = snakemake.params.get("extra", "")

# Setup log
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    """
    (vsearch --{task} {input} {output} {extra}) {log}
    """
)
