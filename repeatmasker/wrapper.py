
from snakemake.shell import shell
import os

print("This is input:", snakemake.input.fa)
print("This is output masked:", snakemake.output.masked)
outdir = os.path.dirname(snakemake.output.masked)
print("This is outdir:", outdir)
extra = snakemake.params.get("extra", "")

compressed_cat = snakemake.output.cat + ".gz"

shell(
    """
    RepeatMasker {extra} -pa {snakemake.threads} {snakemake.input.fa} -dir {outdir}
    if head -n 1 {snakemake.output.out} | grep -q 'There were no repetitive sequences detected'
    then
      ln -sr {snakemake.input.fa} {snakemake.output.masked} && touch {snakemake.output.tbl}
    fi
    if [[ -f {compressed_cat} ]]
    then
      gzip -d {compressed_cat}
    fi
    """
    )