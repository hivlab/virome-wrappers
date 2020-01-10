from snakemake.shell import shell

# Select CD-HIT program, defaults to CD-HIT-EST
program = snakemake.params.get("program", "cd-hit-est")
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Parse threads flag
if snakemake.threads:
    threads = snakemake.threads
else:
    threads = 0

if "psi-cd-est.pl" in program:
    cmd = "({program} -i {snakemake.input} -o {snakemake.output.repres} {extra}) {log}"
else:
    cmd = "({program} -i {snakemake.input} -o {snakemake.output.repres} -T {threads} {extra}) {log}"


shell(cmd)
