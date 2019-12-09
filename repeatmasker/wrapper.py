__author__ = "Taavi Päll"
__copyright__ = "Copyright 2019, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell
import os

outdir = os.path.dirname(snakemake.output[0])
extra = snakemake.params.get("extra", "")

compressed_cat = snakemake.output.cat + ".gz"

shell(
    """
    RepeatMasker {extra} -pa {snakemake.threads} {snakemake.input.fa} -dir {outdir}
    if head -n 1 {snakemake.output.out} | grep -q 'There were no repetitive sequences detected'; then
      ln -sr {snakemake.input.fa} {snakemake.output.masked} && touch {snakemake.output.tbl}
    fi
    if [[ -f {compressed_cat} ]]; then
      gzip -d {compressed_cat}
    fi
    """
)
