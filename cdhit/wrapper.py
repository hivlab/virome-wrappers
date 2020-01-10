from snakemake.shell import shell

# Select CD-HIT program, defaults to CD-HIT-EST
program = snakemake.params.get("program", "cd-hit-est")
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    "({program}"
    " -i {snakemake.input}"
    " -o {snakemake.output.repres}"
    " -T {snakemake.threads}"
    " {extra}) {log}"
)
