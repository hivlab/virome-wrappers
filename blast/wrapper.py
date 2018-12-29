__author__ = "Taavi Päll"
__copyright__ = "Copyright 2018, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbiblastxCommandline

def run_blast(input, output, params):
  """Merge function arguments into dictionary
  """
  options = dict(input)
  options.update(output)
  options.update(params)
  """Compose blastn or blastx command
  """
  if "task" in options:
    cline = NcbiblastnCommandline(**options)
  else:
    cline = NcbiblastxCommandline(**options)
  """Run blast
  """
  stdout, stderr = cline()

if __name__ == '__main__':
    run_blast(snakemake.input, snakemake.output, snakemake.params)
