__author__ = "Taavi Päll"
__copyright__ = "Copyright 2020, Taavi Päll"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"


from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts


java_opts = get_java_opts(snakemake)


shell(
    "picard CreateSequenceDictionary {java_opts} "
    "R={snakemake.input} O={snakemake.output} &> {snakemake.log}"
)
