
from snakemake.shell import shell

reformat_fastq_extra = snakemake.params.get("reformat_fastq_extra", "")
reformat_fasta_extra = snakemake.params.get("reformat_fasta_extra", "")
extra = snakemake.params.get("extra", "")

# Preprocessing command to run.
commands = [
            "reformat.sh in={snakemake.input} out={snakemake.output.fastq} unmappedonly primaryonly {reformat_fastq_extra} {extra}",
            "reformat.sh in={snakemake.output.fastq} out={snakemake.output.fasta} {reformat_fasta_extra} {extra}"
            ]

# Run preprocessing commands.
for cmd in commands:
  shell(cmd)
