import subprocess
from os.path import dirname

params = snakemake.params

taxdict = {}
for k, v in params.items():
    cmd = "get_species_taxids.sh -t {}".format(v)
    process = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, universal_newlines=True)
    taxdict.update({k: process.stdout})

vir, neg = map(lambda keys: {x: taxdict[x] for x in keys}.values(), [["viruses"], list(set(params.keys()) - set(["viruses"]))])

with open(snakemake.output["vir"], "wt") as f:
    f.write(list(vir)[0])

with open(snakemake.output["neg"], "wt") as f:
    f.write("".join(neg))
