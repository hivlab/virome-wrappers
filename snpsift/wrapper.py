__author__ = "Taavi Päll"
__copyright__ = "Copyright 2020, Taavi Päll"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"


from snakemake.shell import shell
extra = snakemake.params.get("extra", "")
fieldnames = snakemake.params.get("fieldnames", "")


shell(
    """
    SnpSift extractFields {extra} {snakemake.input} {fieldnames} > {snakemake.output}
    """
)
